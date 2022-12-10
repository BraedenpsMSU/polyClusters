import re
import main
from time import sleep

from othermain import list_paths
from NCBIQueries import clean_orderdict_to_json
from SearchFactory import NCBIQueryFactory


class QueryTester:
    def __init__(self):
        self.factory = NCBIQueryFactory(config="search_config.json")
        self.delay = .5
        self.pathflag = True
        self.total_size = 0

    def __call__(self, *args, **kwargs):
        current_job = self.factory(*args, **kwargs)
        results_iterator = current_job.get_results()
        sleep(self.delay)
        for code in current_job.run_jobs():
            print(f"{code=}")
            if code == 0:
                _next = next(results_iterator)
                self.total_size += len(str(_next))
                if _next is not None:
                    if not self.pathflag:
                        print(clean_orderdict_to_json(
                            _next['results']
                        ))
                    else:
                        list_paths(_next["results"])
                    print(f"raw_head={str(_next['results'])[:100]}...")
                else:
                    print(_next)
        return current_job.results

help_str ="""
Cmds:
    db1                                 get info on db1
    db1 [<id1>, <id2>, ...]             get info (via fetch) for ids in db1 via fetch
    db1 [<id1>, <id2>, ...] ;summary    get summary of entries with ids in list in db1
    db1 db2                             get a list of links from db1 to db2
    db1 db2 [<id1>, <id2>, ...]         link ids in db1 to entries in db2, does all links
    db1 db2 [<id1>, <id2>, ...] ;<link> link ids in db1 to entries in db2 using specific link
    weak                                toggles error catching
    path                                toggles display mode to list entries as a list of keys
                                            The last value listed is the value return from
                                            the response if you input all values listed before it.
    run                                 Allows you to collect data provided a taxon id
    x                                   exit to python terminal
    h | help                            prints the help menu
"""



if __name__ == '__main__':
    qt = QueryTester()
    results = None
    weakflag = False

    while True:
        try:
            cmd = input('enter command')
            if cmd == 'x':
                break
            if cmd == 'h' or cmd == 'help':
                print(help_str)
                continue
            if cmd == 'weak':
                weakflag = not weakflag
                print('weak flag set - exception will break loop')
                continue
            if cmd == 'path':
                qt.pathflag = not qt.pathflag
                print('path flag set - paths of the json')
                continue
            if cmd == 'run':
                cmd = input('Write a list of Ids to run and commit to the datafiles:')
                if '[' in cmd:
                    ids = cmd.split('[')[1]
                    ids = ids.split(']')[0]
                    ids = list(ids.split(' ')) if ' ' in ids else [ids]
                    main.run(id=ids)
                    continue
            # results = None
            linkname = cmd.split(';')[1].split(' ')[0] if ';' in cmd else None
            # cmd = cmd.replace(',', ' ')
            cmd = cmd.split(';')[0]
            ids = []
            if '[' in cmd:
                ids = cmd.split('[')[1]
                ids = ids.split(']')[0]
                ids = [list(ids.split(' '))] if ' ' in ids else [[ids]]
            # ids = re.findall(r"\[((?:\w|\s)*)]", cmd)
            dbs = re.findall(r"\b[a-zA-z]+\b", cmd)
            print(f"{dbs=}")
            print(f"{ids=}")
            # if ids:
            #     ids = [list(ids[0].split())]
            args = dbs + ids
            print(f"{args=}")
            old_size = qt.total_size
            if linkname:
                if len(dbs) == 2:
                    results = qt(*args, linkname=linkname)
                else:
                    is_summary = linkname == 'summary'
                    results = qt(*args, summary=is_summary)
            else:
                results = qt(*args)
            print(f"Approx_sent_bytes: {qt.total_size:,}")
            print(f"Approx_last_bytes: {qt.total_size - old_size:,}")
        except Exception as e:
            print(e)
            if weakflag:
                raise e
