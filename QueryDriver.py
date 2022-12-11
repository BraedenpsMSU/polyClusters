# Jordan Dood and Braeden Sopp
# 12/5/2022
# This code is to create and run a QueryDriver, a singleton that runs all queries to the DB

# Imports
import time
from SearchFactory import NCBIQueryFactory, NCBICompositeJob
from NCBIQueries import clean_orderdict_to_json
from collections import *


# Static function to check if our result is valid
def validate_result(result):
    """
    :param result: A result from a CompositeJob to be checked \n
    :return: True only if result is of type collections.OrderedDict
    """
    if type(result) == OrderedDict or type(result) is None:
        return True
    else:
        return False


class QueryDriver:
    """
    QueryDriver is a singleton pattern that can be fed inputs, that
    then calls on a factory to return a composite_job object, and
    place that in a que and to run the que use run_query().
    """

# Fields
    _instance = None    # Our one instance of the class
    _job_que = deque()  # Que for storing CompositeJobs
    _job_times = deque()    # A place to store the jobs

# Constants
    WAIT_TIME = 1.1  # Never send more than 3 jobs in 1.1 sec

# Private Methods
    # Overrode __new__() to support singleton pattern
    def __new__(cls, *args, **kwargs):
        if not cls._instance:
            cls._instance = super(QueryDriver, cls).__new__(
                cls, *args, **kwargs)
        return cls._instance

    # Constructor
    def __init__(self):
        try:
            self.factory = NCBIQueryFactory()
        except FileNotFoundError:
            self.factory = NCBIQueryFactory('utils/search_config.json')

    # This is so we can feed jobs into the instance
    def __call__(self, *args, **kwargs):
        self._job_que.append(self.factory(*args, **kwargs))  # TODO Check this call is correct

    # Private method to handle the logic of response to returned codes from a job
    def _return_code_logic(self, code):
        """
        Method to determine the appropriate response for what code we get back
        from running a job.
            :return True if we complete \n
            :return False if not complete, and set a wait on the clock \n
            :raise RuntimeError if we get booted from the server \n
            :raise ValueError if our value is unexpected
        """
        # IF we get None, don't pass anything up and keep going
        if code is None:
            raise ValueError("Code was none")

        # IF we get 0, pass up the data and keep going
        elif code == 0:
            return True

        # IF we get 1, don't pass anything up, wait and keep going
        elif code == 1:
            return False

        # IF we get 2, panic
        elif code == 2:
            raise RuntimeError("CRITICAL FAILURE: kicked off of server!")

        # We should never get anything else, but if we do panic
        else:
            raise ValueError("ERROR: The run_jobs function returned: " + str(code))


    # Private Method to parse returned results
    def _parse_results(self, result):
        """
        :raise ValueError if we were given an invalid result \n
        :return a list of NCBI result objects otherwise
        """
        if validate_result(result):
            result_list = result['results']
            return result_list
        else:
            raise ValueError


   # Private Method to keep us from over pinging the server
    def _proceed(self):
        if len(self._job_times) < 3:
            return
        while time.time() - min(self._job_times) < self.WAIT_TIME:
            pass
        self._job_times.remove(min(self._job_times))
        return



# Public Methods
    # Public Method that yields results to caller as it runs a job and gets them
    def run_query(self):
        """
        This method grabs a CompositeJob and runs it, yielding results if applicable, and ensures no more than 3 pings
        per second to the server, regardless of where it is in the run.
        :yield list of NCBI objects results if possible \n
        :raise GeneratorExit if no jobs to run, ending caller for-loop
        """
        # IF we have a job to run grab it and proceed
        if len(self._job_que) > 0:
            current_job: NCBICompositeJob = self._job_que.popleft()
        # ELSE end the callers for-loop
        else:
            raise GeneratorExit

        # Make sure we start with ok waiting, aka a do for loop thing
        if len(self._job_times) >= 3:
            self._proceed()

        # Run the CompositeJob until done or error
        for return_code, result in current_job.run_jobs():
            self._job_times.append(time.time())
            self._proceed()
            try:
                if self._return_code_logic(return_code):
                    print(current_job.title)
                    try:
                        returnable = self._parse_results(result)
                        yield returnable
                    except ValueError:
                        # We log this error and do nothing right now
                        print("ERROR: A result value in QueryDriver is corrupt!")
                        pass
            except ValueError:
                # We log this error when it is raised, nothing more
                pass
            except RuntimeError:
                # TODO we want to save things here maybe before passing on the error
                raise RuntimeError


# Private methods for testing
def __test_time():
    start = time.time()
    time.sleep(1)
    end = time.time()
    assert end > start


def __test_validation():
    thing = OrderedDict()
    assert validate_result(thing)


def __test_parse_results():
    flag = False
    o_dict = OrderedDict()
    o_dict['results'] = ['thing1', 'thing2']
    test_instance = QueryDriver()
    assert len(test_instance._parse_results(o_dict)) == 2
    try:
        test_instance._parse_results(None)
    except ValueError:
        flag = True
    assert flag


def __test_is_singleton():
    QD1 = QueryDriver()
    QD2 = QueryDriver()
    assert QD1 == QD2


def __test_return_code_logic():
    QD = QueryDriver()
    # Tests for return of 0
    assert QD._return_code_logic(0)
    # Tests for return of 1
    assert not QD._return_code_logic(1)
    # Tests for return of 2
    flag = False
    try:
        QD._return_code_logic(2)
    except RuntimeError:
        flag = True
    assert flag


# Driver for testing
if __name__ == '__main__':
    # __test_time()
    # __test_validation()
    # __test_parse_results()
    # __test_is_singleton()
    # __test_return_code_logic()
    thing1 = QueryDriver()
    thing1('cdd')
    thing1('protein')
    thing1('cdd', 'protein')
    for value in thing1.run_query():
        print(clean_orderdict_to_json(value))
    for value in thing1.run_query():
        print(clean_orderdict_to_json(value))
    for value in thing1.run_query():
        print(clean_orderdict_to_json(value))


