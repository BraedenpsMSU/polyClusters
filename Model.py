from abc import ABC, abstractmethod
from dataclasses import dataclass, asdict

# TODO make Taxon and CDD both inherit a TreeNCBIData class that
#      will allow for easy construction of their Linage
#      models may want to implement a converter method for DBhandler.
#           Do we want this??? Results in some higher coupling.
#           Probably, this is what they are for.
#           Should probably implement a constructor for making instances from DB results


def flatten_with_ids(inp: dict, second=False):
    new_values = dict()
    old_values = []
    for key in inp:
        if isinstance(inp[key], dict):
            if 'db' in inp[key]:
                if not second:
                    inp[key] = inp[key]['uid']
                else:
                    old_values.append(key)
                    new_values[inp[key]['db'] + "_id"] = inp[key]['uid']
            else:
                old_values.append(key)
                new_values.update(flatten_with_ids(inp[key], second=True))
    for key in old_values:
        del inp[key]
    inp.update(new_values)
    return inp

@dataclass(slots=True, frozen=True, order=True)
class NCBIId:
    """
    Stores the unique id of an object in an NCBI database
    """
    db: str
    uid: str


@dataclass(slots=True, frozen=True, eq=True)
class NCBIData(ABC):
    """
    Base class for NCBI data in the project
    """
    id: NCBIId

    def flatten_results(self):
        results_data = asdict(self)
        return flatten_with_ids(results_data)

    def __post_init__(self):
        self.verify()

    @abstractmethod
    def verify(self) -> None:
        """
        Abstract method: insures valid data
        """
        pass


@dataclass(slots=True, frozen=True)
class Domain(NCBIData):
    name: str
    type: str

    def verify(self):
        assert self.id.db == 'cdd'
        assert self.type in ['A', 'B', 'C', 'X', 'Y']
        # TODO assert SuperFamily is never self unless root



@dataclass(slots=True, frozen=True)
class AssemblyGeneData:
    total: int
    protein_coding: int
    non_coding: int
    pseudogene: int

    def gene_count_sum(self):
        return self.protein_coding + self.non_coding + self.pseudogene


@dataclass(slots=True, frozen=True)
class Assembly(NCBIData):
    organism: str
    name: str
    gc_count: int
    length: int
    chromosomes: int
    gene_data: AssemblyGeneData
    taxon_id: NCBIId

    def verify(self):
        assert self.length > self.gc_count > 0
        assert self.chromosomes > 0
        assert self.id.db == 'genome'  # TODO this maybe assembly instead
        assert self.taxon_id.db == 'taxon'
        assert self.gene_data.total >= self.gene_data.gene_count_sum()
        # TODO log where total recorded genes not equal sum


@dataclass(slots=True, frozen=True)
class Taxon(NCBIData):
    super_family_id: NCBIId
    level: str

    def verify(self):
        assert self.id.db == 'taxon'
        assert self.super_family_id.db == 'taxon'
        # TODO assert SuperFamily is never self unless root
        # TODO check that level assigned is appropriate


@dataclass(slots=True, frozen=True)
class Protein(NCBIData):
    name: str
    cdd: Domain

    def verify(self):
        assert self.id.db == 'protein'

# TODO implement gene dataclass
# @dataclass(slots=True)
# class Gene:
#     Uuid: NCBIId
#     name: Optional[str]
#     symbol: Optional[str]
#     size: Optional[int]
#     Taxon: Optional[Taxon]
#     Assemblies: Optional[List[NCBIId]]


if __name__ == "__main__":
    """Demonstrates the verify method above"""
    Id1 = NCBIId(db='cdd', uid='apple')
    Id2 = NCBIId(db='cdd', uid='fruit')
    badId = NCBIId(db='kittens', uid='kittens')
    # domain = Domain(Id1, Id2, 'AppleDomain')
    # domain2 = Domain(Id1, badId, 'badDomain')
    AGD = AssemblyGeneData(3, 1, 1, 1)
    taxon = NCBIId('taxon', '123')
    item = Assembly(NCBIId('genome', '123'),'the_organism', 'test', 22, 44, 1, AGD, taxon)
    test = Domain(NCBIId('cdd', 'thing'),'test','A')
    protein = Protein(NCBIId('protein', 'fdss'), 'name', test)
