import numpy as np
import matplotlib.pyplot as plt
from allensdk.api.queries.rma_api import RmaApi
from allensdk.api.queries.grid_data_api import GridDataApi
import os
from typing import *
import logging


class FISHFetcher:
    def __init__(self) -> None:
        self.rma = RmaApi()
        self.gda = GridDataApi()

    def find_id_ish(self, gene: str, sag_or_cor: str="sagittal",
                    adu_or_dev: str="adult", time_point: str="P56") -> List[int]:
        """Return the ids of Section Data Sets (a single gene experiment)

        Args
        ----
        gene: str
            the gene to search for

        sag_or_cor: str (accepts * wild cards)
            `coronal` or `sagittal` or `*`

        adu_or_dev: str (accepts * wild cards)
            `adult`, `development`, `both`

        time_point: str (it will be autmatically wildcarded)
            e.g. "P56", "E", "E13", "P"

        Returns
        -------
        list of ids:
            sorted by most_recent to mose ancient

        """
        if adu_or_dev == "adult" and "E" in time_point:
            raise ValueError("there is not adult with age %s" % time_point)

        if adu_or_dev == "adult":
            adu_or_dev = "Mouse"
        elif adu_or_dev == "development":
            adu_or_dev = "DevMouse"
        elif adu_or_dev == "both":
            adu_or_dev = "*Mouse"
        else:
            raise ValueError("adu_or_dev='%s' is not valid" % adu_or_dev)
        criteria = ["[failed$eq'false']",
                    "reference_space[name$li'*%s*']" % time_point,
                    "products[abbreviation$li'%s']" % adu_or_dev,
                    "plane_of_section[name$li'%s']" % sag_or_cor,
                    "genes[acronym$eq'%s']" % gene]
        # include='reference_space',
        res = self.rma.model_query("SectionDataSet", criteria=','.join(criteria), only=["id", "qc_date"], num_rows='all')
        if isinstance(res, str):
            raise ValueError("Bad query! Server returned :\n%s" % res)

        qc_date = []
        for i in res:
            if i["qc_date"] is None:
                qc_date.append('')
            else:
                qc_date.append(i["qc_date"])

        ix = np.argsort(qc_date)
        ix = ix[::-1]

        results = []
        for i in ix:
            results.append(int(res[i]["id"]))

        return results

    def download_grid_all(self, gene: str, folder: str='../data', sag_or_cor: str="sagittal",
                          adu_or_dev: str="adult", time_point: str="P56") -> None:
        """Dowloads all the files

         Args
        ----
        gene: str
            the gene to search for

        sag_or_cor: str (accepts * wild cards)
            `coronal` or `sagittal` or `*`

        adu_or_dev: str (accepts * wild cards)
            `adult`, `development`, `both`

        time_point: str (it will be autmatically wildcarded)
            e.g. "P56", "E", "E13", "P"

        """
        ids = self.find_id_ish(gene, sag_or_cor=sag_or_cor, adu_or_dev=adu_or_dev, time_point=time_point) 
        for idd in ids:
            self.gda.download_expression_grid_data(idd,
                                                   path=os.path.join(folder,
                                                                     "%s_%s_%s_%s.zip" % (gene,sag_or_cor,time_point,idd)))

    def download_grid_recent(self, gene: str, folder: str='../data', sag_or_cor: str="sagittal",
                             adu_or_dev: str="adult", time_point: str="P56") -> None:
        """Dowloads the most recently qc-ed file among the ones available

         Args
        ----
        gene: str
            the gene to search for

        sag_or_cor: str (accepts * wild cards)
            `coronal` or `sagittal` or `*`

        adu_or_dev: str (accepts * wild cards)
            `adult`, `development`, `both`

        time_point: str (it will be autmatically wildcarded)
            e.g. "P56", "E", "E13", "P"

        """
        ids = self.find_id_ish(gene, sag_or_cor=sag_or_cor, adu_or_dev=adu_or_dev, time_point=time_point)
        try:
            idd = ids[0]
            output_path = os.path.join(folder, "%s_%s_%s_%s.zip" % (gene, sag_or_cor, time_point, idd))
            self.gda.download_expression_grid_data(idd, path=output_path)
        except IndexError:
            logging.warn("Experiment %s was never performed" % gene)
            pass
