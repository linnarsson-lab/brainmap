import numpy as np
import matplotlib.pyplot as plt
from allensdk.api.queries.rma_api import RmaApi
from allensdk.api.queries.grid_data_api import GridDataApi
import os
import glob
from typing import *
import brainmap as bm
import re
import logging
from collections import OrderedDict


class LimitedSizeDict(OrderedDict):
    def __init__(self, *args: Any, **kwds: Any) -> None:
        self.size_limit = kwds.pop("size_limit", None)  # type: int
        OrderedDict.__init__(self, *args, **kwds)
        self._check_size_limit()

    def __setitem__(self, key: Any, value: Any) -> None:
        OrderedDict.__setitem__(self, key, value)
        self._check_size_limit()

    def _check_size_limit(self) -> None:
        if self.size_limit is not None:
            while len(self) > self.size_limit:
                self.popitem(last=False)


class ISHFetcher:
    ''' A downloader object for Section Data Sets

    Methods
    -------

    find_id_ish:
        Returns the ids of Section Data Sets (a single gene experiment) sorted by qc time

    download_grid_all:
        Dowloads all the expression energy 3d density file (200um grid) that satisfy the query

    download_grid_recent:
        Dowloads the most recently qc-ed expression energy 3d density file (200um grid) that satisfy the query

    Attributes
    ----------
    rma:
        Rma Api instance
    gda
        GridData Api instance
    res
        results of the find_id_ish query
    '''
    def __init__(self) -> None:
        self.rma = RmaApi()
        self.gda = GridDataApi()
        self.res = None  # type: List

    def find_id_ish(self, gene: str, sag_or_cor: str="sagittal",
                    adu_or_dev: str="adult", time_point: str="P56") -> List:
        """Returns the ids of Section Data Sets (a single gene experiment)

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
        self.res = self.rma.model_query("SectionDataSet", criteria=','.join(criteria), only=["id", "qc_date"], num_rows='all')
        if isinstance(self.res, str):
            raise ValueError("Bad query! Server returned :\n%s" % self.res)

        if self.res == []:
            return []

        qc_date = []
        for i in self.res:
            if i["qc_date"] is None:
                qc_date.append('')
            else:
                qc_date.append(i["qc_date"])

        ix = np.argsort(qc_date)
        ix = ix[::-1]

        results = []
        for i in ix:
            results.append(int(self.res[i]["id"]))

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
                                                                     "%s_%s_%s_%s.zip" % (gene, sag_or_cor, time_point, idd)))

    def download_grid_recent(self, gene: str, folder: str='../data', sag_or_cor: str="sagittal",
                             adu_or_dev: str="adult", time_point: str="P56") -> Union[str, bool]:
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

        Returns
        -------
        output_path: output_path or bool
            if the download was successfull returns the path to the file otherwise False

        """
        ids = self.find_id_ish(gene, sag_or_cor=sag_or_cor, adu_or_dev=adu_or_dev, time_point=time_point)
        try:
            idd = ids[0]
            output_path = os.path.join(folder, "%s_%s_%s_%s.zip" % (gene, sag_or_cor, time_point, idd))
            self.gda.download_expression_grid_data(idd, path=output_path)
            return output_path
        except IndexError:
            logging.warn("Experiment %s was never performed" % gene)
            return False


class ISHLoader:
    def __init__(self, root: str, adu_or_dev: str="adult",
                 time_point: str="P56", priority: List[str]=["coronal", "sagittal"]) -> None:
        self.root = root
        self.adu_or_dev = adu_or_dev
        self.time_point = time_point
        self.priority = priority
        self.regex_list = [re.compile("(.*?)%s(.*?).zip" % i) for i in self.priority]
        assert os.path.isdir(self.root), "%s is not a folder" % self.root
        self._fetcher = ISHFetcher()
        self.index = {}  # type: Dict[str, str]
        self._build_index()
        self._cache = LimitedSizeDict(size_limit=300)  # type: LimitedSizeDict

    def _build_index(self) -> None:
        """Build dict gene->file
        """
        self.duplicates = []  # type: List[str]
        file_list = glob.glob(os.path.join(self.root, "*_*.zip"))
        if file_list == []:
            return
        all_genes, all_paths = zip(*[(os.path.split(i)[-1].split("_")[0], i) for i in file_list])
        for n, regex in enumerate(self.regex_list[::-1]):
            for i, (path, gene) in enumerate(zip(all_paths, all_genes)):
                if regex.match(path):
                    if gene in self:
                        self.duplicates.append(gene)
                    self.index[all_genes[i]] = path
        logging.debug("%i duplicates were found" % len(self.duplicates))

    def __contains__(self, value: Any) -> bool:
        return value in self.index

    def __getitem__(self, value: str) -> np.ndarray:
        if value in self._cache:
            return self._cache[value]
        elif value in self:
            path = self.index[value]
            vol_data = bm.AllenVolumetricData(filename=path)
            self._cache[value] = vol_data
            return vol_data
        else:
            logging.debug("%s was not in root, attempting dowload" % value)
            for sag_or_cor in self.priority:
                output_path = self._fetcher.download_grid_recent(gene=value, folder=self.root, sag_or_cor=sag_or_cor,
                                                                 adu_or_dev=self.adu_or_dev, time_point=self.time_point)
                if output_path and isinstance(output_path, str):
                    logging.debug("%s slicing derived dataset was found" % sag_or_cor)
                    self.index[value] = output_path
                    return bm.AllenVolumetricData(filename=self.index[value])
            raise KeyError("gene %s is not available in root or for dowload in the Allen Brain Atlas" % value)
