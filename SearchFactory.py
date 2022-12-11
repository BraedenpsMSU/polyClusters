import collections
import json
from collections import OrderedDict
from functools import lru_cache, partial
from typing import Iterable, Sequence, Optional, TextIO, List, Final

import NCBIQueries
from datetime import datetime, timezone

import NCBIJobPacket
from Model import NCBIId
from NCBIJobPacket import NCBIBaseJobPacket

NCBI_CONFIG: Final[str] = "search_config.json"


class NCBICompositeJob:
    """
    Collection NCBI Queries that must be run to get the finish some task.

    :ivar title: Name of collection. Should indicate its goal of task
    :ivar run_date: Date at which the job was run. Set to None before instance is run.
    :ivar results: Collection of results from Jobs. Set to None before instance is run.
    :ivar finished: True if job has been run to completion successfully else False.
    """

    def __init__(self, title: str, jobs: list[NCBIBaseJobPacket]):
        """
        :param title: title of this job collection.
        :param jobs: list of Queries that Inherit BaseNCBIQuery to be run.
        """
        self.title: str = title
        self._jobs: Sequence[NCBIBaseJobPacket] = tuple(jobs)
        self._complete: int = 0
        self.finished: bool = False
        self.run_date: Optional[datetime] = None
        self.last_request_time: Optional[datetime] = None
        self.results_yielded: int = 0
        self.results_produced: int = 0
        self._results: list[OrderedDict] = []
        self.results: Optional[list[OrderedDict]] = None

    def __len__(self):
        return len(self._jobs)

    def run_jobs(self) -> Iterable[int]:
        """
        Produces generator that yields:
            0 if result was successful, \n
            1 if result was unsuccessful due to server being down \n
            2 if result was unsuccessful due to exceed API ratelimit \n
        Only proceeds to next job if success.

        :return: A int corresponding to if computing the results was successful.
        :raises: RuntimeError when a request is sent that cause the server to HTTP status other
                 than 200, 438, and 503
        """
        self._results = list()
        self.run_date = datetime.now(timezone.utc)
        self.results_yielded: int = 0
        self.results_produced: int = 0
        self._complete: int = 0
        self._results: list[OrderedDict] = []
        self.finished: bool = False
        while self._complete < len(self):
            runner = self._jobs[self._complete]
            for query in runner.run():
                self.last_request_time = runner.get_last_request_time()
                if not query.success \
                        and not query.exceed_limit \
                        and not query.server_unavailable:
                    if query.status_code == 500:
                            raise LookupError("Server is having internal issues")
                    raise RuntimeError(f"Bad Server request with https-code {query.status_code}")
                elif not query.success \
                        and query.exceed_limit:
                    yield 2, None
                elif not query.success:
                    yield 1, None
                else:
                    runner.progress()
                    current_result = None
                    if runner.is_finished():
                        current_result = self.format_result(
                            runner.generate_results(),
                            time_of_request=self.last_request_time
                        )
                        assert current_result['results'] is not None
                        self._results.append(
                            current_result
                        )
                    yield 0, current_result
            self._complete += 1
        self.results = OrderedDict()
        for result in self._results:
            self.results.update(result['results'])
        self.results = self.format_result(
            self.results,
            time_of_request=self.last_request_time
        )
        self.finished = True

    def format_result(self, result, time_of_request=None):
        outp = OrderedDict()
        outp['composite_star_time'] = self.run_date
        outp['time_of_last_request'] = time_of_request
        outp['command'] = self.title
        outp['results'] = result
        return outp


class NCBIQueryFactory:

    def __init__(self, config=NCBI_CONFIG):
        self.config_file: str = config
        f: TextIO
        with open(self.config_file, 'r') as f:
            self.config: dict = json.load(f)

    def __call__(self, *args, **kwargs) -> NCBICompositeJob:
        """
        Jank but fast implementation \n
        0 Argument Call:
            self() gets general entrez info \n
        1 Argument Calls
            self(arg0: str) gets info db passed in \n
            self(arg0: Sequence[NCBIId]) gets information on all ids given \n
        2 Arguments Calls
            self(arg0: str, arg1: str) gets link info db args0 to db args1 \n
            self(arg0: str, arg1: Sequence[str]) get info on all ids in arg1 in database arg0 \n
            self(arg0: str, arg1: Sequence[NCBIId]) gets all link data from an NCBIId \n
        3 Argument Calls
            self(arg0: str, arg1: str, arg2: Sequence[str])
                arg0: db you are linking from \n
                arg1: db you are linking to \n
                arg2: list of id for the link operation \n
        """
        if kwargs.get('taxon'):
            if kwargs.get("page_token"):
                return self.taxon_query(kwargs.get('taxon'), page_token=kwargs.get("page_token"))
            return self.taxon_query(kwargs.get('taxon'))
        if kwargs.get('lineage'):
            return self.lineage_query(kwargs.get('lineage'))
        if len(args) == 0:
            return self.getinfo(args)
        if len(args) == 1:
            if isinstance(args[0], str):
                return self.getinfo(*args)
            elif isinstance(args[0], Sequence):
                return self.get_data(*args, **kwargs)
            else:
                raise ValueError("Bad Arguments")
        if len(args) == 2:
            assert isinstance(args[0], str), "First argument must be string when 2 arguments are given"
            search_obj = args[1]
            if isinstance(search_obj, str):
                return self.get_links_info(*args)
            elif isinstance(search_obj, Sequence):
                if isinstance(search_obj[0], str):
                    new_input: list[NCBIId] = []
                    for _id in search_obj:
                        new_input.append(NCBIId(args[0], _id))
                    if kwargs.get('summary'):
                        return self.get_Summary(new_input)
                    return self.get_data(new_input)
                return self.get_link(*args, **kwargs)
        if len(args) == 3:
            assert isinstance(args[0], str), "First argument must be string when 3 arguments are given"
            assert isinstance(args[1], str), "Second argument must be string when 3 arguments are given"
            assert isinstance(args[2], Sequence), "Third argument must be a Sequence when 3 arguments are given"
            new_input: list[NCBIId] = []
            for _id in args[2]:
                new_input.append(NCBIId(args[0], _id))
            return self(args[1], new_input, **kwargs)
        else:
            raise ValueError('To many arguments!')

    def getinfo(self, *param) -> NCBICompositeJob:
        if param[0] and param[0] not in self.config['db']:
            raise ValueError(f"{repr(param[0])} is not a valid DB")
        _job = NCBIQueries.EntrezInfo() if not param[0] else NCBIQueries.EntrezInfo(db=param[0])
        value = NCBIJobPacket.NCBIMonoJobPacket(
            _job,
            extractor=lambda x: x["eInfoResult"]
        )
        _title = "Info" if not param[0] else f"Info {param[0]}"
        return NCBICompositeJob(title=_title, jobs=[value])

    def get_links_info(self, *args):
        def linkname_extractor(inp: OrderedDict, filter_term=args[0]) -> OrderedDict:
            outp = OrderedDict()
            for item in inp["eInfoResult"]["DbInfo"]['LinkList']['Link']:
                if item['DbTo'] == filter_term:
                    outp[item['Name']] = item
            return outp
        assert len(args) == 2, "To many args have been given"
        if args[0] not in self.config['db'] or args[1] not in self.config['db']:
            raise ValueError(f"{args[0]} or {args[1]} is not a valid DB")
        _job = NCBIQueries.EntrezInfo(db=args[0])
        value = NCBIJobPacket.NCBIMonoJobPacket(
            _job,
            extractor=partial(linkname_extractor, filter_term=args[1])
        )
        _title = f"LinkInfo {args[0]} {args[1]}"
        return NCBICompositeJob(title=_title, jobs=[value])

    def get_link(self, param: str, param1: Sequence[NCBIId], linkname: Optional[str] = None):
        if len(param1) < self.config['chunk']:
            _jobs = [NCBIJobPacket.NCBIMonoJobPacket(
                    NCBIQueries.EntrezLink(param1[0].db, param, [ncbi_id.uid for ncbi_id in param1], linkname=linkname)
                )
            ]
            _title = f"link {param1[0].db} {param} {[ncbi_id.uid for ncbi_id in param1]}"
        else:
            _jobs = [NCBIJobPacket.NCBIMonoJobPacket(
                NCBIQueries.EntrezLink(param1[0].db, param, [ncbi_id.uid for ncbi_id in chunk], linkname=linkname)
                ) for chunk in self.list_chunker(param1)
            ]
            _title = f"link {param1[0].db} {param} {[ncbi_id.uid for ncbi_id in param1]}"
        return NCBICompositeJob(_title, _jobs)

    def get_data(self, *param, **kwargs):
        if len(param[0]) <= self.config["chunk"]:
            _jobs = [NCBIJobPacket.NCBIMonoJobPacket(
                    NCBIQueries.EntrezFetch(db=param[0][0].db, id=[ncbi_id.uid for ncbi_id in param[0]],
                                            html_tags=self.config["tags"][param[0][0].db])
                )
            ]
            _title = f"Fetch Data from {param[0][0].db} {[ncbi_id.uid for ncbi_id in param[0]]}"
        else:
            _jobs = [NCBIJobPacket.NCBIMonoJobPacket(
                    NCBIQueries.EntrezFetch(db=param[0][0].db, id=[ncbi_id.uid for ncbi_id in chunk],
                                            html_tags=self.config["tags"][param[0][0].db])
                ) for chunk in self.list_chunker(param[0])
            ]
            _title = f"Fetch Data from {param[0][0].db} chunking {[ncbi_id.uid for ncbi_id in param[0]]}"
        return NCBICompositeJob(_title, _jobs)

    def get_Summary(self, *param, **kwargs):
        if len(param[0]) <= self.config["chunk"]:
            _jobs = [NCBIJobPacket.NCBIMonoJobPacket(
                    NCBIQueries.EntrezSummary(db=param[0][0].db, id=[ncbi_id.uid for ncbi_id in param[0]])
                )
            ]
            _title = f"Summary from {param[0][0].db} chunking {[ncbi_id.uid for ncbi_id in param[0]]}"

        else:
            _jobs = [NCBIJobPacket.NCBIMonoJobPacket(
                NCBIQueries.EntrezSummary(db=param[0][0].db, id=[ncbi_id.uid for ncbi_id in chunk])
                ) for chunk in self.list_chunker(param[0])
            ]
            _title = f"Summary from {param[0][0].db} chunking {[ncbi_id.uid for ncbi_id in param[0]]}"
        return NCBICompositeJob(_title, _jobs)


    def taxon_query(self, taxon, page_token=None):
        if page_token:
            _jobs = [NCBIJobPacket.NCBIMonoJobPacket(
                NCBIQueries.GeneDatasetQuery(taxon, page_token=page_token)
            )
            ]
        else:
            _jobs = [NCBIJobPacket.NCBIMonoJobPacket(
                NCBIQueries.GeneDatasetQuery(taxon)
            )
            ]
        _title = f"getting genes"
        return NCBICompositeJob(_title, _jobs)

    def lineage_query(self, param):
        _jobs = [NCBIJobPacket.NCBIMonoJobPacket(
            NCBIQueries.TaxonSummary(param)
        )
        ]
        _title = f"getting taxons"
        return NCBICompositeJob(_title, _jobs)

    def list_chunker(self, input_list: Sequence[NCBIId], chunk_size=None) -> List[List[NCBIId]]:
        """
        Returns a list of list of NCBI such that each list in the output is
        no larger than the chunk_size. \n
        If chunk_size is not set then we have chunk_size is set to the chunk value of the config file. \n
        Base on the grouper method in the itertools recipes. \n
        See: https://docs.python.org/3/library/itertools.html#itertools-recipes \n


        :param input_list:
        :param chunk_size:
        :return: List of List of NCBI such that each list less than chunk_size.
        """
        chunk_size = self.config['chunk'] if chunk_size is None else chunk_size
        _iter = [iter(input_list)] * chunk_size
        output = []
        while True:
            pack = [next(_iter_item, None) for _iter_item in _iter]
            pack = list(filter(lambda x: x is not None, pack))
            if pack:
                output.append(pack)
            else:
                break
        return output


if __name__ == "__main__":
    pass
    # factory = NCBIQueryFactory()
    # test1 = factory()
    # test2 = factory('cdd')
    # for i in test1.run_jobs():
    #     print(i)
    # for i in test1.get_results():
    #     print(i)
    #     results1 = i
    # for j in test2.run_jobs():
    #     print(j)
    # for j in test2.get_results():
    #     print(j)
    #     results2 = j
    # small_input_list = (
    #     NCBIId('cdd', '1'),
    #     NCBIId('cdd', '2')
    # )
    # test3 = factory('cdd', 'protein')
    # for j in test3.run_jobs():
    #     print(j)
    # for j in test3.get_results():
    #     print(j)
    #     results3 = j
    # print(NCBIQueries.clean_orderdict_to_json(results3['results']))
    # # large_input_list = tuple([str(i) for i in range(50)])
    # # large_ncbiid_input_lisst = tuple([NCBIId('cdd', i) for i in large_input_list])
    # # test3 = factory(
    # #     'cdd',
    # #     large_input_list
    # # )
    f = NCBIQueryFactory()
    ptest = ['167637656', '167637389', '167634873', '167631881', '167539906',
             '167532179', '167529623', '167514896', '167514205', '167511999',
             '167454746', '167405094', '167385332', '167384251', '167383931',
             '167363943', '167344430', '167218155', '167214966', '167211248',
             '167069206', '167058337', '167015347', '167013667', '167007523',
             '166996884', '166993941', '166993558', '166985422', '166984312',
             '166982709', '166031757', '166028210', '165933451', '165908538',
             '165900780', '165900578', '165899947', '165894112', '165870887',
             '165868953', '165868883', '164715679', '164715609', '164713396',
             '164687503', '164665154', '164603277', '163942629', '163940679',
             '163939037', '163864826', '163862876', '163861234', '163842025',
             '163814308', '163812765', '163798979', '163790960', '163789371',
             '163762723', '163759711', '163746181', '163741508', '163738780',
             '163691876', '162955501', '162454878', '162452637', '162447983',
             '162436623', '162283308', '161986090', '161388108', '161385242',
             '161380926', '161165460', '161163219', '160947116', '160944134',
             '160940361', '160933951', '160893795', '160888657', '160886580',
             '160885861', '160878978', '160427644', '160356656', '159875225',
             '159873802', '159043840', '159037713', '158972986', '158449943',
             '158446250', '158444810', '158436941', '158338749', '158315803',
             '158308990']
    _ptest2 = ['167637656', '167637389', '167634873', '167631881', '167539906',
              '167532179', '167529623', '167514896', '167514205', '167511999',
              ]
    id_values = [NCBIId('protein', _id) for _id in _ptest2]
    chunks = f.list_chunker(id_values, 5)
    test_comp = NCBICompositeJob('title',
                     [
                         NCBIJobPacket.NCBIMonoJobPacket(
                             NCBIQueries.EntrezSummary(
                                 db=chunk[0].db,
                                 id=[item.uid for item in chunk]
                             )
                         ) for chunk in chunks
                     ])


