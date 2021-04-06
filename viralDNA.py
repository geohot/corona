import copy

class viralDNA:
    """
    Python provides its own interface of Prototype via `copy.copy` and
    `copy.deepcopy` functions. And any class that wants to implement custom
    implementations have to override `__copy__` and `__deepcopy__` member
    functions.
    """

    # Defaults are set for corona virus
    def __init__(self, untranslated_region_f = 265, orf1a_i = 266,orf1a_f = 13483, orf1b_i = 13468,orf1b_f = 21555, spike_gp_i = 21563, spike_gp_f = 25384, orf3a_i = 25393, orf3a_f = 26220, envelope_p_i = 26245, envelope_p_f = 26472, membrane_gp_i = 26523, membrane_gp_f = 27191, orf6_i = 27202, orf6_f = 27387, orf7a_i = 27394, orf7a_f = 27759, orf7b_i = 27756, orf7b_f = 27887, orf8_i = 27894, orf8_f = 28259, n_p_i = 28274, n_p_f = 29533, orf10_i = 29558, orf10_f = 29674):
        # in front "the untranslated leader sequence that ends with the Transcription Regulation Sequence"
        self.untranslated_region = cc[0 : untranslated_region_f]

        self.orf1a = translate(cc[orf1a_i - 1 : orf1a_f], True)

        # cc[266-1+4398*3:13468] = 'TTT_TTA_AAC' aka 'X_XXY_YYZ'
        # https://en.wikipedia.org/wiki/Ribosomal_frameshift
        # Programmed −1 Ribosomal Frameshifting
        # TODO: add this to the translate function with automatic detection
        self.orf1b = translate(cc[orf1b_i - 1 : orf1b_f], False).strip("*")  # chop off the stop, note this doesn't have a start

        # exploit vector, this attaches to ACE2. also called "surface glycoprotein"
        # https://www.ncbi.nlm.nih.gov/Structure/pdb/6VYB -- open state
        # https://www.ncbi.nlm.nih.gov/Structure/pdb/6VXX -- closed state
        # 1273 amino acids
        #   S1  = 14-685
        #   S2  = 686-1273
        #   S2' = 816-1273
        # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2750777/
        self.spike_glycoprotein = translate(cc[spike_gp_i - 1 :spike_gp_f], True)

        # Forms homotetrameric potassium sensitive ion channels (viroporin) and may modulate virus release.
        self.orf3a =  translate(cc[orf3a_i-1:orf3a_f], True)

        # These two characteristics stick out, used in assembly aka they package the virus
        self.envelope_protein = translate(cc[envelope_p_i - 1 : envelope_p_f], True)  # also known as small membrane
        self.membrane_glycoprotein = translate(cc[membrane_gp_i - 1 : membrane_gp_f], True)
        self.orf6 = translate(cc[orf6_i - 1 : orf6_f], True)
        self.orf7a = translate(cc[orf7a_i - 1 : orf7a_f], True)
        self.orf7b = translate(cc[orf7b_i - 1 : orf7b_f], True)  # is this one real?
        self.orf8 = translate(cc[orf8_i - 1 : orf8_f], True)

        # https://en.wikipedia.org/wiki/Capsid
        # Packages the positive strand viral genome RNA into a helical ribonucleocapsid
        # Includes the "internal" protein (from Coronavirus Pathogenesis)
        # https://www.sciencedirect.com/topics/veterinary-science-and-veterinary-medicine/human-coronavirus-oc43
        self.nucleocapsid_phosphoprotein = translate(cc[n_p_i - 1 : n_p_f], True)

        # might be called the internal protein (Coronavirus Pathogenesis)
        self.orf10 = translate(cc[orf10_1 - 1 : orf10_f], True)

    def __copy__(self):
        """
        Shallow copy. 
        Call with copy.copy 
        Return value of copy.copy is new shallow copy
        """

        # Copies of the class fields (nested objects)
        untranslated_region = copy.copy(self.untranslated_region)
        orf1a = copy.copy(self.orf1a)
        orf1b = copy.copy(self.orf1b)
        spike_glycoprotein = copy.copy(self.spike_glycoprotein)
        orf3a = copy.copy(self.orf3a)
        envelope_protein = copy.copy(self.envelope_protein)
        membrane_glycoprotein = copy.copy(self.membrane_glycoprotein)
        orf6 = copy.copy(self.orf6)
        orf7a = copy.copy(self.orf7a)
        orf7b = copy.copy(self.orf7b)
        orf8 = copy.copy(self.orf8)
        nucleocapsid_phosphoprotein = copy.copy(self.nucleocapsid_phosphoprotein)
        orf10 = copy.copy(self.orf10)

        new = self.__class__(
            self.some_int, some_list_of_objects, some_circular_ref
        )
        new.__dict__.update(self.__dict__)

        return new

    def __deepcopy__(self, memo={}):
        """
        Deep copy
        Call with copy.deepcopy 
        Return value of copy.deepcopy is new deep copy
        Memo is the dictionary used by the `deepcopy` library - prevents circular references.
        """
        # Copies of the class fields (nested objects)
        untranslated_region = copy.deepcopy(self.untranslated_region, memo)
        orf1a = copy.deepcopy(self.orf1a, memo)
        orf1b = copy.deepcopy(self.orf1b, memo)
        spike_glycoprotein = copy.deepcopy(self.spike_glycoprotein, memo)
        orf3a = copy.deepcopy(self.orf3a, memo)
        envelope_protein = copy.deepcopy(self.envelope_protein, memo)
        membrane_glycoprotein = copy.deepcopy(self.membrane_glycoprotein, memo)
        orf6 = copy.deepcopy(self.orf6, memo)
        orf7a = copy.deepcopy(self.orf7a, memo)
        orf7b = copy.deepcopy(self.orf7b, memo)
        orf8 = copy.deepcopy(self.orf8, memo)
        nucleocapsid_phosphoprotein = copy.deepcopy(self.nucleocapsid_phosphoprotein, memo)
        orf10 = copy.deepcopy(self.orf10, memo)

        # Then, let's clone the object itself, using the prepared clones of the
        # nested objects.
        new = self.__class__(
            self.untranslated_region, orf1a, orf1b, spike_glycoprotein, orf3a, envelope_protein, membrane_glycoprotein, orf6, orf7a, orf7b, orf8, nucleocapsid_phosphoprotein, orf10
        )
        new.__dict__ = copy.deepcopy(self.__dict__, memo)

        return new


if __name__ == "__main__":
    custom_viralDNA = viralDNA()
    shallow_copied_viralDNA = copy.copy(custom_viralDNA)
    deep_copied_viralDNA = copy.deepcopy(custom_viralDNA)

    