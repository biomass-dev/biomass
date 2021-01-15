import os
import re
import numpy as np
from dataclasses import dataclass
from typing import List, Tuple

from .template import Model


class BioMassModel(Model):
    pass


@dataclass
class ExecModel(object):
    model: BioMassModel

    def get_individual(self, paramset: int) -> np.ndarray:
        best_generation = np.load(
            os.path.join(
                self.model.path,
                "out",
                f"{paramset:d}",
                "generation.npy",
            )
        )
        best_individual = np.load(
            os.path.join(
                self.model.path,
                "out",
                f"{paramset:d}",
                f"fit_param{int(best_generation):d}.npy",
            )
        )
        return best_individual

    def load_param(self, paramset: int) -> Tuple[list, list]:
        best_individual = self.get_individual(paramset)
        (x, y0) = self.model.sp.update(best_individual)
        return x, y0

    def get_executable(self) -> List[int]:
        n_file = []
        try:
            fitparam_files = os.listdir(
                os.path.join(
                    self.model.path,
                    "out",
                )
            )
            for file in fitparam_files:
                if re.match(r"\d", file):
                    n_file.append(int(file))
            empty_folder = []
            for i, nth_paramset in enumerate(n_file):
                if not os.path.isfile(
                    os.path.join(
                        self.model.path,
                        "out",
                        f"{nth_paramset:d}",
                        "generation.npy",
                    )
                ):
                    empty_folder.append(i)
            for i in sorted(empty_folder, reverse=True):
                n_file.pop(i)
        except FileNotFoundError as e:
            print(e)
        return n_file