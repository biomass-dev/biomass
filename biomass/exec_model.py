class ExecModel(object):
    def __init__(self, model):
        self.model_path = model.__path__[0]
        self.parameters = model.C.NAMES
        self.species = model.V.NAMES
        self.rxn = model.ReactionNetwork()

    def get_properties(self):
        """ Get number of model reactions, species and parameters
        """
        biological_processes = self.rxn.group()
        model_reactions = 0
        for reactions_in_process in biological_processes:
            model_reactions += len(reactions_in_process)
        model_species = len(self.species)
        model_parameters = len(self.parameters)

        print(
            "{:d} reactions\n{:d} species\n{:d} parameters".format(
                model_reactions, model_species, model_parameters
            )
        )