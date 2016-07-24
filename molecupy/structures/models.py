from .molecules import AtomicStructure

class Model(AtomicStructure):

    def __init__(self):
        pass


    def __getattr__(self, attribute):
        if attribute == "_atoms":
            return set()
        else:
            return self.__getattribute__(attribute)
