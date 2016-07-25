class Pdb:

    def __init__(self, data_file):
        self._data_file = data_file


    def data_file(self):
        return self._data_file


    def classification(self):
    	return self._data_file.classification()


    def deposition_date(self):
    	return self._data_file.deposition_date()


    def pdb_code(self):
    	return self._data_file.pdb_code()


    def is_obsolete(self):
    	return self._data_file.is_obsolete()


    def obsolete_date(self):
    	return self._data_file.obsolete_date()


    def replacement_code(self):
    	return self._data_file.replacement_code()


    def title(self):
    	return self._data_file.title()


    def split_codes(self):
    	return self._data_file.split_codes()


    def caveat(self):
    	return self._data_file.caveat()


    def keywords(self):
    	return self._data_file.keywords()


    def experimental_techniques(self):
    	return self._data_file.experimental_techniques()


    def model_count(self):
    	return self._data_file.model_count()


    def model_annotations(self):
    	return self._data_file.model_annotations()


    def authors(self):
    	return self._data_file.authors()


    def revisions(self):
    	return self._data_file.revisions()


    def supercedes(self):
    	return self._data_file.supercedes()


    def supercede_date(self):
    	return self._data_file.supercede_date()


    def journal(self):
    	return self._data_file.journal()
