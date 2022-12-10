import json
from abc import ABC, abstractmethod
from collections import OrderedDict
from typing import Callable, Optional, Final

import xmltodict as x2d
from requests import get, Response


def clean_orderdict_to_json(inp: OrderedDict) -> str:
    return json.dumps(inp, indent=4)


class BaseNCBIQuery(ABC):
    """
    Base class for Entries object.

    Will have 2 main immediate children
        1) Entrez Search
        2) NCBI DataSet REST API search
    """

    def __init__(self, *args, **kwargs):
        self.status_code: Optional[int] = 0
        self.cmd: Callable[[], Response] = get
        self.success: Optional[bool] = None
        self.exceed_limit: Optional[bool] = None
        self.server_unavailable: Optional[bool] = None

    @staticmethod
    def kwargs_to_html(kwargs, start=False):
        if not kwargs:
            return ''
        args_list = [f"{key}={kwargs[key]}" for key in kwargs]
        lead = "?" if start else "&"
        return lead + '&'.join(args_list)

    @abstractmethod
    def _get_base_api_url(self) -> str:
        """
        Provide string that represents API server url

        :return: Server URL
        """
        pass

    @abstractmethod
    def _get_cmd_string(self) -> str:
        """
        Provides of command part of REST cmd.

        :return: String of the command
        """
        pass

    @abstractmethod
    def _process_response(self, response: Response) -> OrderedDict:
        """
        Takes the processed request and give JSON version on response back.

        :return:
        """
        pass

    @abstractmethod
    def _determine_success(self, response: Response):
        """
        determine if the response from the server was successful.
        Set instance variables: success and exceed_limit and server_unavailable
        exceed_limit or server_unavailable implies not success.
        server_unavailable and exceed_limit cannot both be true
        """
        pass

    def request(self) -> OrderedDict:
        """
        Make request to server.
        :return: OrderDictionary of the JSON response of search
        """
        query = self._get_base_api_url() + self._get_cmd_string()  # creates query
        response = self.cmd(query)  # runs query
        self.status_code = response.status_code
        self._determine_success(response)  # checks if query was successful
        if self.success:  # if successful
            return self._process_response(response)  # returns json as OrderDict of response
        return OrderedDict()  # else return empty OrderDict


# noinspection SpellCheckingInspection
class BaseEntrezQuery(BaseNCBIQuery, ABC):
    def __init__(self, *args, **kwargs):
        super(BaseEntrezQuery, self).__init__(*args, **kwargs)
        self.cmd = get
        self.html_tags = self.kwargs_to_html(kwargs['html_tags']) if kwargs.get('html_tags') is not None else ''

    def _get_base_api_url(self) -> str:
        return "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

    def _determine_success(self, response: Response):
        if response.status_code != 200:  # check if OKAY
            self.success = False
            # if not OKAY check if limit exceeded.
            self.exceed_limit = True if response.status_code == 429 else False
            # then check if the server is unavailable as well
            self.server_unavailable = True if response.status_code == 503 else False
        else:
            self.success = True
            self.exceed_limit = False
            self.server_unavailable = False


class EntrezInfo(BaseEntrezQuery):
    entrez_cmd: Final[str] = 'einfo'

    def __init__(self, *args, **kwargs):
        super(EntrezInfo, self).__init__(*args, **kwargs)
        self.db = kwargs.get('db')

    def _get_cmd_string(self) -> str:
        db_param = '' if self.db is None else f"?db={self.db}"
        return f"{self.__class__.entrez_cmd}.fcgi{db_param}"

    def _process_response(self, response: Response) -> OrderedDict:
        return x2d.parse(response.content.decode(response.apparent_encoding))


class EntrezFetch(BaseEntrezQuery):
    entrez_cmd: Final[str] = 'efetch'

    def __init__(self, *args, **kwargs):
        super(EntrezFetch, self).__init__(*args, **kwargs)
        self.db = kwargs.get('db')
        self.ids = kwargs.get('id')  # TODO can be None!
        if not self._valid_input():
            raise ValueError

    def _get_cmd_string(self) -> str:
        db_param = f"?db={self.db}"
        ids_param = f"&id={','.join(self.ids)}"
        return f"{self.__class__.entrez_cmd}.fcgi{db_param}{ids_param}{self.html_tags}"

    def _process_response(self, response: Response) -> OrderedDict:
        return x2d.parse(response.content.decode(response.apparent_encoding))

    # TODO: can this return json? see, want xml v2 instead?
    #  https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly

    def _valid_input(self):
        flag = True
        if self.db is None:  # Check we have a db
            print("EtrezFetch does not have required db parameter")
            flag = False
        if len(self.ids) >= 200:  # Check we don't ask for too much
            print("EtrezFetch is too long")
            flag = False
        if type(self.ids) != list or len(self.ids) < 1:  # Check we have a list of ids
            print("EtrezFetch was not given the required list of IDs")
            flag = False
        return flag
    # TODO check are the ids str or obj? Ever done without list?


class EntrezSummary(BaseEntrezQuery):  # Jordan working here!!
    entrez_cmd: Final[str] = 'esummary'

    def __init__(self, *args, **kwargs):
        super(EntrezSummary, self).__init__(*args, **kwargs)
        self.db = kwargs.get('db')
        self.ids = kwargs.get('id')  # TODO check this is the input name, is it ever none?
        if not self._valid_input():
            raise ValueError

    def _get_cmd_string(self) -> str:
        db_param = f"?db={self.db}"
        ids_param = f"&id={','.join(self.ids)}"
        return f"{self.__class__.entrez_cmd}.fcgi{db_param}{ids_param}&retmode=json{self.html_tags}"

    def _process_response(self, response: Response) -> OrderedDict:
        return response.json(object_pairs_hook=OrderedDict)

    # https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly

    def _valid_input(self):
        flag = True
        if self.db is None:  # Check we have a db
            print("EtrezSummary does not have required db parameter")
            flag = False
        if len(self.ids) >= 200:  # Check we don't ask for too much
            print("EtrezSummary is too long")
            flag = False
        if not isinstance(self.ids, list) or len(self.ids) < 1:  # Check we have a list of ids
            print("EtrezSummary was not given the required list of IDs")
            flag = False
        return flag
    # TODO check are the ids str or obj? Ever done without list?


class EntrezLink(BaseEntrezQuery):
    entrez_cmd: Final[str] = 'elink'

    def __init__(self, *args, **kwargs):
        self.source_db: str = args[0]
        self.foreign_db: str = args[1]
        self.id: list[str] = args[2]
        self.linkname: str = kwargs.get("linkname")
        super(EntrezLink, self).__init__(args[3:])

    def _get_cmd_string(self) -> str:
        id_arg = ','.join(self.id)
        link_param = '' if self.linkname is None else f"&linkname={self.linkname}"
        return f"{self.__class__.entrez_cmd}.fcgi?dbfrom={self.source_db}&db={self.foreign_db}&id={id_arg}{link_param}&retmode=json"

    def _process_response(self, response: Response) -> OrderedDict:
        #return x2d.parse(response.content.decode(response.apparent_encoding))
        return response.json(object_pairs_hook=OrderedDict)


class EntrezPost(BaseEntrezQuery):
    entrez_cmd: Final[str] = 'epost'

    def __init__(self, *args, **kwargs):
        super(EntrezPost, self).__init__(*args, **kwargs)
        self.entrez_cmd: str = args[0]
        self.db: str = args[1]
        self.ids: list[str] = args[2]
        self.kwargs = kwargs

    def _get_cmd_string(self) -> str:
        return f"{self.__class__.entrez_cmd}.fcgi?db={self.db}&id={','.join(self.ids)}"

    def _process_response(self, response: Response) -> OrderedDict:
        results = x2d.parse(response.content.decode(response.apparent_encoding))
        return results['ePostResult']


class EntrezReceivePost(BaseEntrezQuery):
    def __init__(self, *args, **kwargs):
        self.entrez_cmd: str = args[0]
        self.db = args[1]
        self.query_id: str = args[2]['QueryKey']
        self.web_env: str = args[2]['WebEnv']
        self.kwargs = self.kwargs_to_html(kwargs, start=False)
        # WARNING jank stuff ahead. Determining how to process the response at runtime
        if self.entrez_cmd == EntrezLink.entrez_cmd:
            setattr(self, "_process_response", getattr(EntrezLink, "_process_response"))
        elif self.entrez_cmd == EntrezSummary.entrez_cmd:
            setattr(self, "_process_response", getattr(EntrezSummary, "_process_response"))
        elif self.entrez_cmd == EntrezFetch.entrez_cmd:
            setattr(self, "_process_response", getattr(EntrezFetch, "_process_response"))
        else:
            raise ValueError(f"{self.entrez_cmd} is not a valid Entrez Command for a post operation")
        super(EntrezReceivePost, self).__init__(args[3:], **dict())

    def _get_cmd_string(self) -> str:
        return f"{self.entrez_cmd}.fcgi?db={self.db}&query_key={self.query_id}&WebEnv={self.web_env}{self.kwargs}"

    def _process_response(self, response: Response) -> OrderedDict:
        raise NotImplementedError("This should have been bound on creation!")


class GeneDatasetQuery(BaseNCBIQuery):
    def __init__(self, *args, **kwargs):
        super(GeneDatasetQuery, self).__init__(*args, **kwargs)
        self.taxon_id = args[0]
        self.tags = '&page_size=999'
        if kwargs.get('page_token'):
            self.tags += f"&page_token={kwargs.get('page_token')}"

    def _get_cmd_string(self) -> str:
        return f"/gene/taxon/{self.taxon_id}?table_fields=gene-id{self.tags}"

    def _process_response(self, response: Response) -> OrderedDict:
        return response.json(object_pairs_hook=OrderedDict)

    def _determine_success(self, response: Response):
        if response.status_code != 200:  # check if OKAY
            self.success = False
            # if not OKAY check if limit exceeded.
            self.exceed_limit = True if response.status_code == 429 else False
            # then check if the server is unavailable as well
            self.server_unavailable = True if response.status_code == 503 else False
        else:
            self.success = True
            self.exceed_limit = False
            self.server_unavailable = False

    def _get_base_api_url(self) -> str:
        return "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"

    # TODO

# https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/208964"

class TaxonSummary(BaseNCBIQuery):
    # TODO
    def __init__(self, *args, **kwargs):
        super(TaxonSummary, self).__init__(*args, **kwargs)
        self.taxon_id = args[0]

    def _get_cmd_string(self) -> str:
        return f"/taxonomy/taxon/{','.join(self.taxon_id)}"

    def _process_response(self, response: Response) -> OrderedDict:
        return response.json(object_pairs_hook=OrderedDict)

    def _determine_success(self, response: Response):
        if response.status_code != 200:  # check if OKAY
            self.success = False
            # if not OKAY check if limit exceeded.
            self.exceed_limit = True if response.status_code == 429 else False
            # then check if the server is unavailable as well
            self.server_unavailable = True if response.status_code == 503 else False
        else:
            self.success = True
            self.exceed_limit = False
            self.server_unavailable = False

    def _get_base_api_url(self) -> str:
        return "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"


class DataSetAssemblySummary(BaseNCBIQuery):
    # TODO
    pass


if __name__ == "__main__":
    x = EntrezInfo()
