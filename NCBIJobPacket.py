import collections
from abc import ABC, abstractmethod
from collections import OrderedDict
from typing import Optional, OrderedDict, Tuple, Iterable, Callable

from NCBIQueries import BaseNCBIQuery
from datetime import datetime, timezone


class NCBIBaseJobPacket(ABC):
    def __init__(self, *args, **kwargs):
        self._finished: bool = False
        self._last_ping_timestamp: Optional[datetime] = None

    @abstractmethod
    def run(self) -> Iterable[BaseNCBIQuery]:
        pass

    @abstractmethod
    def generate_results(self) -> OrderedDict:
        pass

    @abstractmethod
    def progress(self):
        pass

    @abstractmethod
    def __len__(self):
        pass

    def is_finished(self) -> bool:
        return self._finished

    def get_last_request_time(self) -> datetime:
        assert self._last_ping_timestamp is not None, "A job must be run to make a request"
        return self._last_ping_timestamp


class NCBIMonoJobPacket(NCBIBaseJobPacket):
    def __init__(self, *args, **kwargs):
        self.__query: BaseNCBIQuery = args[0]
        self.__result_extractor: Callable[[OrderedDict], OrderedDict] = kwargs.get('extractor', (lambda x: x))
        self.__result: Optional[OrderedDict] = None
        super(NCBIMonoJobPacket, self).__init__()

    def run(self) -> Iterable[BaseNCBIQuery]:
        while not self._finished:
            self._last_ping_timestamp = datetime.now(timezone.utc)
            self.__result = self.__query.request()
            yield self.__query

    def progress(self):
        self._finished = True

    def generate_results(self) -> OrderedDict:
        assert self.__result is not None, \
            "Job must finish running before results are available"
        return self.__result_extractor(self.__result)

    def __len__(self):
        return 1


class NCBIFollowUpJobPacket(NCBIBaseJobPacket):
    def __init__(self, *args, **kwargs):
        self.__query: BaseNCBIQuery = args[0]
        self.__follow_up_query: Optional[BaseNCBIQuery] = None
        self.__follow_up: Callable[[BaseNCBIQuery, OrderedDict], BaseNCBIQuery] = args[1]
        self.__result_extractor: Callable[[OrderedDict, OrderedDict], OrderedDict] = kwargs.get(
            'extractor', (lambda r1, r2: r2)
        )
        self.__result: Optional[OrderedDict] = None
        self.__follow_up_result: Optional[OrderedDict] = None
        self.state = 1
        super(NCBIFollowUpJobPacket, self).__init__()

    def run(self) -> Iterable[BaseNCBIQuery]:
        while not self._finished:
            if self.state == 1:
                self._last_ping_timestamp = datetime.now(timezone.utc)
                self.__result = self.__query.request()
                yield self.__query
            if self.state == 2:
                assert self.__follow_up_query is not None, "Must progress the job to run!"
                self._last_ping_timestamp = datetime.now(timezone.utc)
                self.__follow_up_result = self.__follow_up_query.request()
                yield self.__follow_up_query

    def generate_results(self) -> OrderedDict:
        assert self.__result is not None and self.__follow_up_result is not None, \
            "Job must finish running before results are available"
        return self.__result_extractor(self.__result, self.__follow_up_result)

    def progress(self):
        if self.state == 2:
            self._finished = True
        if self.state == 1:
            assert self.__result is not None, "Must run job to progress!"
            self.__follow_up_query = self.__follow_up(self.__query, self.__result)
            self.state += 1

    def __len__(self):
        return 2


class NCBICollectorJobPacket(NCBIBaseJobPacket):
    def __init__(self, *args: NCBIBaseJobPacket, **kwargs):
        assert len(args) >= 2, "At least two packets are required for a collector"
        assert all(isinstance(arg, NCBIBaseJobPacket) for arg in args), "Every argument must be a packet"
        self.__queries: Tuple[NCBIBaseJobPacket] = args
        self.__accum_results: OrderedDict = collections.OrderedDict()
        self.__size = sum(len(query) for query in self.__queries)
        self.__current_job_index = 0
        super(NCBICollectorJobPacket, self).__init__(*args, **kwargs)

    def run(self) -> Iterable[BaseNCBIQuery]:
        while self.__current_job_index < len(self.__queries):
            for query in self.__queries[self.__current_job_index]:
                self._last_ping_timestamp = self.__queries[self.__current_job_index].get_last_request_time()
                yield query

    def generate_results(self) -> OrderedDict:
        assert self.is_finished(), "Job must be finished be results are generated"
        return self.__accum_results

    def progress(self):
        if not self._finished:
            if self.__queries[self.__current_job_index].is_finished():
                self.__accum_results.update(
                    self.__queries[self.__current_job_index].generate_results()
                )
                self.__current_job_index += 1
                if self.__current_job_index == len(self.__queries):
                    self._finished = True
            else:
                self.__queries[self.__current_job_index].progress()

    def __len__(self):
        return self.__size
