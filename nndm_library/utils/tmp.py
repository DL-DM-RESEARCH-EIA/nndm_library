
class ReadLhe(ReadFileBase):
    def __init__(self, path, variables_interest=["energy", "ps"], catch=1, verbose=1):
        """
        ReadLhe: Class
        -------------------------------------------------------------
        * Parameters:
        - variables_interest: List. Variable to extract from the lhe.
        egs. ["energy","angle"], ["energy","px","py"] .... By default: ["energy","angle"]
            KeyArgs:
            - path: String. Path to search for the txt with the path for lhe's
            - catch: Int. Identificate the particle. 1 is the outgoing electron (measured) . Default 1
        * Output
        - DataFrame with all readed events.
        """
        ReadFileBase.__init__(self, path) 

        self.doc = ET.parse(path)
        self.root = self.doc.getroot()

        # catch is used to select the electrons that go away
        self.catch = catch
        self.cartessian = ["px", "py", "pz"]

        self.dic_var = {}
        if "ps" in variables_interest:
            variables_interest = variables_interest + self.cartessian
        for var in [v for v in variables_interest if v != "ps"]:
            self.dic_var[var] = np.array([])

        self.verbose = verbose
        self._read()

    def _getAngle(self, momentums, axis=2):
        """
        Calculate the angle of the particles starting from a list of the form [px, py, pz].
        This with respect to the "axis" element.
        """
        momentumTotal = np.linalg.norm(momentums)
        angule = np.arccos(momentums[axis] / momentumTotal)
        return angule

    def _read(self):
        """
        procced with certain LHE reading
        """

        if self.verbose:
            decide_tqdm_apply = tqdm.tqdm
        else:
            decide_tqdm_apply = lambda x: x

        for bloque in decide_tqdm_apply(self.root):
            Evento = []
            if bloque.tag == "event":  # Evaluate if data comes from valid event
                datos = bloque.text
                for x in datos.splitlines():
                    a = x.split()
                    if len(a) == 13:
                        Evento.append(a)
                energia_electrones = [float(Evento[1][9]), float(Evento[3][9])]
                momentum_electrones = [
                    [float(Evento[1][i]) for i in [6, 7, 8]],
                    [float(Evento[3][i]) for i in [6, 7, 8]],
                ]  # [[px,py,pz]]
                for var in self.dic_var.keys():
                    if var in self.cartessian:  # Append all momentums
                        # self.dic_var[var] = np.append(self.dic_var[var],momentum_electrones[self.catch])
                        dif = self.cartessian.index(var)
                        self.dic_var[self.cartessian[dif]] = np.append(
                            self.dic_var[self.cartessian[dif]],
                            momentum_electrones[self.catch][dif],
                        )
                    elif "angle" == var:  # Append the angle
                        self.dic_var[var] = np.append(
                            self.dic_var[var],
                            self._getAngle(momentum_electrones[self.catch], axis=2),
                        )
                    elif "energy" == var:  # Append the energy :: 'energy'
                        self.dic_var[var] = np.append(
                            self.dic_var[var], energia_electrones[self.catch]
                        )
                    elif "P" == var:  # Get magnitude of Four Momentum :: P
                        self.dic_var[var] = np.append(
                            self.dic_var[var],
                            np.linalg.norm(momentum_electrones[self.catch]),
                        )

        self.dic_var = pd.DataFrame.from_dict(self.dic_var)

    def get(self):
        return self.dic_var

    def get_momenta(self):
        """
        name: name of the file, which is assumed to contain the values of the eps and the mass to be extracted from it
        format example is as follows: eta_decay_events_mk_0.38_eps2_5.404557191441203e-07.lhe

        return: momentum and energy
            px, py, pz, energy
        """
        return self.dic_var["energy"], self.dic_var["px"], self.dic_var["py"], self.dic_var["pz"]