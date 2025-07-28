#!/usr/bin/env python

"""
Produce a version of GTDB compatible with MG-TK

This will attempt to download given version of GTDB, and format files to the 
required structure for MG-TK.
There are some differences in how GTDB provide files between version,
so the standard function is written for 220, and future or past version
might require some custom work. However I have tried to write this in an
extensible manner, so you only need to subclass and override a few methods
where custom processing is required, and the commandline interface can be
reused.

If adding a new version of GTDB, then you will need to do a few things
depending on what has changed:
* URLs for files have changed - override GTDBVersion.__init__ to provide
  the new URLs in 'required_urls'. See example GTDB111 class.
* Organisation of marker genes in tars has changed - override _extract_markers
  to handle different storage structure
* Organisation of metadata has changed - override _format taxonomy.
In each of the methods, downloaded_file is passed in, giving you a dictionary
of the downloaded files you can access, so you will be able to get to any
file among those downloaded in each method. See subclass GTDB111 for example.

You will always needs to add to the GTDBVersion.from_version() method to
determine when to use your new subclass.

If adding to this script, please try to use only modules from the standard
library to make it easier to run on different systems. Also please use type
annotations at least in function definitions (arguments, return types).
"""

# TODO: Ask Jogi if BenchPro can handle 220

__author__ = "Anthony Duncan"
__version__ = 0.1

# Define settings for different GTDB version
import argparse
from collections import defaultdict, namedtuple
import csv
import gzip
import io
import json
import logging
import os
import sys
import shutil
from statistics import mode
import subprocess
from datetime import datetime
import tarfile
from typing import Callable, Dict, List, NamedTuple, Optional, Set, Tuple, Literal, Any
from pathlib import Path

# Configure logging
logger: logging.Logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
ch: logging.Handler = logging.StreamHandler()
fm: logging.Formatter = logging.Formatter(
    fmt='%(levelname)s [%(asctime)s]: %(message)s',
    datefmt='%d/%m/%Y %I:%M:%S'
)
ch.setFormatter(fm)
logger.addHandler(ch)

class GTDBVersion:
    """
    Class to represent a version of GTDB and methods to download and
    format it for use in MG-TK. All default implementations are intended
    to work for v220. For versions which need alternative implementation,
    subclass and override methods.

    Most methods assume you are running on a UNIX like systems with normal
    utilities (gunzip, mv, zcat, etc) available.
    """

    GTDBTK_VERSIONS: Dict[float, Tuple[str, str]] = {
        226.0: ("2.4.0", "Current"),
        220.0: ("2.4.0", "Current"),
        214.0: ("2.1.0", "2.3.2"),
        207.2: ("2.1.0", "2.3.2"),
        207.0: ("2.0.0", "2.0.0"),
        202.0: ("1.5.0", "1.7.0")
    }

    def __init__(
            self,
            version: float,
            base_url: str = ("https://data.ace.uq.edu.au/public/gtdb/data"
                             "/releases/"),
            required_urls: Optional[Dict[str, str]] = None,
            test: bool = False
    ) -> None:
        """
        :param version: GTDB version
        :param base_url: Base URL for GTDB download. Provided incase GTDB
        moves some day
        :param required_urls: Files which need to be download. Provide as
        a dictionary with the key being a label and the value the url.
        Either full URLs or suffixes to base_url. If this is left blank,
        files correct at v220 will be added to this list. See __default_urls
        for details of what the standard URLs are.
        :param test: Make dummy files instead of downloasing, to allow testing
        script logic without downloading large files.
        """
        self.version: float = version
        self.base_url: str = base_url
        self.required_urls: Dict[str, str] = self.__default_urls(required_urls)
        self.__dl_name: Literal["wget", "curl", "test"]
        self.__dl_cmd: Callable
        self.__set_dl_cmd(test)
    
    def __set_dl_cmd(self, test: bool = False) -> None:
        """Determine which utility to use for file downloads. Prefer wget.
        :param test: Make dummy files instead of download."""
        if test:
            self.__dl_name = "test"
            self.__dl_cmd = cmd_test
            logger.warning(
                "Using test mode for downloads. This will create dummy files "
                "rather than performing downloads."
            )
        elif self.program_avail("wget"):
            self.__dl_name = "wget"
            self.__dl_cmd = cmd_wget
        elif self.program_avail("curl"):
            self.__dl_name = "curl"
            self.__dl_cmd = cmd_curl
            logger.warning(
                "Using curl for downloads. This does not support directory" \
                "downloads so will not be able to download GTDBtk database" \
                "in parts."
            )
        else:
            logger.error("wget or curl must be available to perform downloads")
            raise Exception("Missing wget or curl")

    def __default_urls(
            self,
            required_urls: Optional[Dict[str, str]]
            ) -> Dict[str, str]:
        """
        Make standard files to download assuming v220 format.

        Default keys are 'bac_markers', 'bac_metadata', 'bac_taxonomy',
        'arc_markers', 'arc_metadata', 'arc_taxonomy', and 'tk_data'.
        Any labels prefixed 'tk_' should be GTDBtk specific data, and will
        be skipped by default (they are ~110GB).
        """
        if required_urls is None:
            dir: str = (f"{self.base_url}/release{self.version:.0f}"
                        f"/{self.version:.1f}")
            return dict(
                bac_markers=(
                    f"{dir}/genomic_files_all/bac120_marker_genes_all_r"
                    f"{self.version:.0f}"
                    ".tar.gz"
                ),
                bac_taxonomy=(
                    f"{dir}/bac120_taxonomy_r{self.version:.0f}.tsv.gz"),
                bac_metadata=(
                    f"{dir}/bac120_metadata_r{self.version:.0f}.tsv.gz"),
                arc_markers=(
                    f"{dir}/genomic_files_all/ar53_marker_genes_all_r"
                    f"{self.version:.0f}"
                    ".tar.gz"
                ),
                arc_taxonomy=(
                    f"{dir}/ar53_taxonomy_r{self.version:.0f}.tsv.gz"),
                arc_metadata=(
                    f"{dir}/ar53_metadata_r{self.version:.0f}.tsv.gz"),
                tk_database=(
                    f"{dir}/auxillary_files/gtdbtk_package/full_package/"
                    f"gtdbtk_r{self.version:.0f}_data.tar.gz"
                ),
                tk_database_parts=(
                    f"{dir}/auxillary_files/gtdbtk_package/split_package/")
            )
        return required_urls

    @staticmethod
    def program_avail(cmd: str) -> bool:
        """Determine if a given program is available on shell."""
        subp: subprocess.CompletedProcess = subprocess.run(
                cmd, 
                shell=True,
                capture_output=True
            )
        if subp.returncode == 0:
            return True
        if "--help" in str(subp.stderr):
            return True
        return False

    def download(
            self,
            download_dir: Path,
            error_existing: bool = False,
            offline: bool = False,
            download_tk: Optional[Literal[False, "full", "split"]] = False
        ) -> Dict[str, Path]:
        """
        Download required files to a download directory. File will be retained
        after download. Downloads are done by calling wget.

        :param download_dir: Directory to save files to
        :param error_existing: Whether to error if local files exist. If this
        is False (default), it will be assumed the local files are correct
        and will not be replaced.
        :param offline: Do not attempt to download data. Will not check if files
        are complete.
        :param download_tk: Download GTDBtk data. This is skipped by default
        due to size.
        :return: Mapping of url to local location
        """
        dl_check: Path = download_dir / ".download.finished"
        if dl_check.exists():
            logger.info(
                "Download completion marker found (%s). If you want to redo "
                "downloads delete this file.", dl_check)
            downloaded = {
                label: download_dir / (
                    label if is_url_dir(url)
                    else url.split("/")[-1]
                ) 
                for label, url in self.required_urls.items()
            }
            return downloaded

        logger.info("Downloads required:")
        # Remove GTDBtk files if not requested
        tk_remove: List[str] = ["tk_database", "tk_database_parts"]
        if download_tk == "full":
            tk_remove = ["tk_database_parts"]
        elif download_tk == "split":
            tk_remove = ["tk_database"]

        filt_urls: Dict[str, str] = {
            label: url for label, url in self.required_urls.items()
            if label not in tk_remove
        }
        for label, url in filt_urls.items():
            logger.info("* %s\t%s", label, url)
        downloaded: Dict[str, Path] = {}

        for label, url in filt_urls.items():
            is_dir: bool = is_url_dir(url)
            fname: str = label if is_dir else url.split("/")[-1]
            local: Path = download_dir / fname
            cmd: str = self.__dl_cmd(url, local)

            # Determine if file exists
            if local.exists() and error_existing:
                raise FileExistsError()
            
            # If in offline mode, add to downloaded if exists
            if offline:
                if not local.exists():
                    logger.error("Extract mode, but %s not found", str(local))
                    raise FileNotFoundError(local)
                logger.info("Extract mode, download skipped for %s", label)
                downloaded[label] = local
                continue

            logger.info("Download %s", label)
            logger.debug(
                "Download %s\nSource: %s\nDestination: %s\nCommand: %s",
                label,
                url,
                str(local),
                cmd
            )

            # Not offline, so attempt download
            subp: subprocess.CompletedProcess = subprocess.run(
                cmd, 
                shell=True,
                capture_output=True
            )
            if subp.returncode != 0 or not local.exists():
                # Check if we got a n error due to fully received
                if not "fully retrieved" in subp.stderr.decode('utf-8'):
                    logger.error(
                        "Download of %s failed. Dumping output.", url)
                    logger.error("Stderr")
                    logger.error(subp.stderr)
                    logger.error("Stdout")
                    logger.error(subp.stdout)
                    self.__clean_failed_gtdbtk_dir(local)
                    raise Exception(subp)
            downloaded[label] = local
        if not offline:
            self._processing_metadata(
                dest=download_dir,
                dir=str(download_dir),
                tk=str(download_tk)
            )
            dl_check.write_text(str(datetime.now()))
        return downloaded

    def __clean_failed_gtdbtk_dir(
            self,
            dir: Path
    ) -> None:
        """Remove any incomplete parts of the GTDBtk data from the download
        directory. wget does not resume these properly, so we want to only
        keep parts we can ensure are complete."""

        logger.info(
            "Directory download failed. Removing incomplete file "
            "as cannot resume directory download."
        )
        gzs: List[Path] = list(dir.glob("*.tar.gz*"))
        mode_size: float = mode(gz.stat().st_size for gz in gzs)
        last_file: Path = max(gzs, key=lambda x: x.stat().st_mtime)
        logger.debug("Last file: %s", last_file)
        logger.debug("Mode size: %s", str(mode_size))
        if (
            (last_file.stat().st_size != mode_size)
            or
            (len(gzs) < 3)
        ):
            # The last downloaded file is too small, or we don't have enough
            # parts to know what size the parts should be, so delete and redo
            # the last one
            logger.info("Deleting %s", last_file)
            last_file.unlink()

    def __extract_seqs(
            self,
            src: Path,
            dest: Path
        ) -> List[Path]:
        """
        Extract all FAA/FNA from a tar to a destination, optionally ignoring
        files which already exist. Returns a list of the extracted FNA/FAA
        files. This is split out into a method as archaea and bacteria need
        to be extracted separately, and file for some shared marker genes
        concatenated.

        :param dest: Extraction directory
        :returns: List of all file extracted
        """
        logger.debug("Source %s", src)
        logger.debug("Destination: %s", dest)
        dest.mkdir(exist_ok=True, parents=True)
        
        # Memory usage and speed of tarfile is poor seemingly, so directly
        # using tar in shell

        # --strip-components flattens the file structure (remove leading 2
        # parts of the path)
        cmd: str = (
            f"tar -xf {str(src)} --strip-components 2 -C {str(dest)}"
        )
        logger.debug("Comand: %s", cmd)

        subp: subprocess.CompletedProcess = subprocess.run(
            cmd, 
            shell=True,
            capture_output=True
        )
        subp.check_returncode()
        logger.debug("Extraction complete")
        return list(dest.glob("*.fna")) + list(dest.glob("*.faa"))
    
    def _extract_markers(
            self,
            dest: Path,
            downloaded_files: Dict[str, Path]
    ) -> None:
        """
        Extract marker genes to the MG-TK format. Not intended to be called
        directly, but for future GTDB version with different structure,
        it may be useful to be able to override this method. Default extract
        method opens 'labels' in download files and put faa and fna in the
        same directory.

        :param dest: Destination directory. Contents will be overwritten.
        :param downloaded_files: Files downloaded for parsing, with key being
        a label and value path of local file.
        """

        # .marker.finished file inidicates this step has alredy completed
        finished_marker: Path = dest / ".marker.finished"
        logger.debug("Check for finished marker file: %s", finished_marker)
        if finished_marker.exists():
            logger.info("Marker extraction already complete.")
            return
    
        logger.info("Extracting archaeal markers")
        arc_dest: Path = dest / "arc"
        arc_dest.mkdir(parents=True, exist_ok=True)
        arc_markers: List[Path] = self.__extract_seqs(
            src=downloaded_files['arc_markers'],
            dest=arc_dest,
        )
        logger.info("Extracting bacterial markers")
        bac_dest: Path = dest
        bac_markers: List[Path] = self.__extract_seqs(
            src=downloaded_files['bac_markers'],
            dest=bac_dest
        )

        logger.info("Combining archaeal and bacterial markers")
        # Determine which markers exist for both archaea and bacteria
        # Where there are multiple, concatenate them
        marker_files: Dict[str, list] = defaultdict(list)
        for file in bac_markers + arc_markers:
            marker_files[file.name].append(file)
        for f_name, files in (
            (n, f) for n, f in marker_files.items() if len(f) > 1):
            bac_f, arc_f = files[0], files[1]
            cmd: str = f"cat {arc_f} >> {bac_f} && rm {arc_f}"
            logger.debug(f"Concatenate marker {f_name}")
            logger.debug(f"Command: {cmd}")
            proc_res: subprocess.CompletedProcess = subprocess.run(
                cmd,
                shell=cmd,
                capture_output=True
            )
            proc_res.check_returncode()
        # Move any remaining archaeal marker files to the main directory
        if len(list(arc_dest.glob("*.f*"))) > 1:
            logger.info("Moving remaining archaea specific markers")
            cmd: str = f"mv {arc_dest}/*.f* {bac_dest}"
            logger.debug(f"Command: {cmd}")
            proc_res: subprocess.CompletedProcess = subprocess.run(
                cmd,
                shell=True,
                capture_output=True
            )
            if proc_res.returncode != 0:
                logger.error(proc_res.stderr)
                logger.error(proc_res.stdout)
            proc_res.check_returncode()
        arc_dest.rmdir()
        finished_marker.write_text(str(datetime.now()))

    def _extract_tk(
        self,
        dest: Path,
        downloaded_files: Dict[str, Path]
    ) -> None:
        """
        Extract GTDBtk data

        :param dest: Destination for extracted files
        :param download_files: Source files downloaded
        """

        finished_marker: Path = dest / ".tk.finished"
        if finished_marker.exists():
            logger.info("GTDBtk extracation already completed.")
            return

        parts: bool = "tk_database_parts" in downloaded_files.keys()
        src: Path = downloaded_files['tk_database']
        dest.mkdir(parents=True, exist_ok=True)
        cmd: str = f"tar -xf {src} -C {dest}"
        if parts:
            logger.info("GTDBtk database parts will be concatenated for "
            "extraction. Concatenated version will be deleted after extraction "
            "and parts retained")
            part_src: Path = downloaded_files['tk_database_parts']

            # Check if concatenation completed
            concat: Path = src
            src_size: float = sum(
                p.stat().st_size for p in part_src.glob("*gz*part*"))
            dest_size: float = (
                0 if not concat.exists() else concat.stat().st_size
            )
            cat_cmd: str = (
                f"cat {part_src / '*.tar.gz.part*'} > "
                f"{concat}"
            )
            # Allow a tolerance in filesize sum
            TOL = 0.0001
            if not (
                (dest_size > (src_size * (1 - TOL)))
                and
                (dest_size < (src_size * (1 + TOL)))
            ):
                logger.info("Concatenating GTDBtk database parts to %s", concat)
                logger.debug("Destination: %s", concat)
                logger.debug("Command: %s", cat_cmd)

                proc_res: subprocess.CompletedProcess = subprocess.run(
                    cat_cmd,
                    shell=True,
                    capture_output=True
                )
                proc_res.check_returncode()

                logger.info("GTDBtk concatenation complete")
                
            cmd = (
                f"tar -xf {concat} -C {dest} && "
                f"rm -f {concat}"
            )
        
        logger.info("Extracting GTDBtk database to %s", dest)
        logger.debug("Source: %s", src)
        logger.debug("Destination: %s", dest)
        logger.debug("Command: %s", cmd)

        proc_res: subprocess.CompletedProcess = subprocess.run(
            cmd,
            shell=True,
            capture_output=True
        )
        proc_res.check_returncode()

        logger.info("GTDBtk extraction complete")
        finished_marker.write_text(str(datetime.now()))

    def _format_taxonomy(
            self,
            dest: Path,
            downloaded_files: Dict[str, Path],
            taxonomy_headers: List[str] = ["accession", "taxonomy"]
    ) -> None:
        """
        Make lineage and clustering tables in tab separated format. Not
        intended to be called directly, can be overriden for files which
        are not in the v220 format. This uses standard library functions to
        make it easier to use as part of generic install process.
        """

        logger.info("Producing taxonomy tables")
        logger.debug("Destination: %s", dest)

        logger.info("Concatenate source archaeal and bacterial tables")
        bac_tax, arc_tax = (
            downloaded_files['bac_taxonomy'], downloaded_files['arc_taxonomy'])
        bac_md, arc_md = (
            downloaded_files['bac_metadata'], downloaded_files['arc_metadata'])
        comb_md: Path = bac_md.parent / "combined_md.tsv.gz"
        comb_tax: Path = bac_md.parent / "combined_tax.tsv.gz"
        # Can't get process substitution working via subprocess, so having
        # to make some intermediate files (and don't want to handle in this 
        # script opening to append etc)
        cmd: str = (
            f"cat {bac_tax} {arc_tax} > {comb_tax} && "
            f"zcat {arc_md} | tail -n +2 | gzip > tmp_arc_md.tsv.gz && "
            f"cat {bac_md} tmp_arc_md.tsv.gz > {comb_md} && "
            "rm tmp_arc_md.tsv.gz"
        )
        logger.debug(f"Command: {cmd}")
        proc_res: subprocess.CompletedProcess = subprocess.run(
            cmd,
            capture_output=True,
            shell=True
        )
        proc_res.check_returncode()

        logger.info("Making species lineage")
        species: Set[str] = set()
        with (
            gzip.open(comb_tax, 'rt') as tbl,
            open(
                dest / f"gtdb_r{self.version:.0f}_lineageGTDB.tab",
                "w") as lin_o
        ):
            lineage_write = csv.writer(
                lin_o, delimiter="\t", quoting=csv.QUOTE_NONE)
            for row in csv.DictReader(
                tbl, delimiter="\t", fieldnames=taxonomy_headers):
                accession, lineage = row['accession'], row['taxonomy']
                lineage_parts: List[str] = lineage.split(";")
                # Only add a species once
                if lineage_parts[-1] in species:
                    continue
                species.add(lineage_parts[-1])
                lineage_write.writerow(
                    [lineage_parts[-1], 0, *lineage_parts]
                )
        
        logger.info("Making species to representative accession mapping")
        with (
            gzip.open(comb_md, 'rt') as tbl,
            open(
                dest / f"gtdb_r{self.version:.0f}_clustering.tab",
                "w") as clust_o
        ):
            lineage_write = csv.writer(
                clust_o, delimiter="\t", quoting=csv.QUOTE_NONE)
            for row in csv.DictReader(tbl, delimiter="\t"):
                if row['gtdb_representative'] == "t":
                    tax: str = row['gtdb_taxonomy']
                    accession: str = row['accession']
                    species_str = tax.split(";")[-1]
                    lineage_write.writerow([species_str, accession])

        comb_md.unlink()
        comb_tax.unlink()

        logger.info("Taxonomy tables finished")

    def _accessions_from_markers(self, marker_dir: Path) -> Set[str]:
        """
        Find all accession which relate to extracted marker sequences.

        In r207, all marker genes related to representative species. This seems
        to not be the case in subsequent releases. This method will find all
        accession which related to any marker sequences, so their taxonomy
        can be extracted.

        :param marker_dir: Path of extracted marker gene sequences
        """
        
        # Parse FAA as shorter and same contents
        accs: Set[str] = set()
        for faa in marker_dir.glob("**/*.faa"):
            logger.debug("Get accessions for %s", faa)
            for line in (x for x in faa.open("rt") if x[0] == ">"):
                acc: str = line[1:].strip()
                accs.add(acc)
        return accs

    @staticmethod
    def ask_config_update(do_update: Optional[bool] = None) -> bool:
        """Determine whether or not to do MG-TK configuration update."""

        if do_update is not None:
            return do_update

        if not sys.__stdin__.isatty():
            logger.info(
                "Run in non-interactive terminal, skipping config update")
            return False

        mgtk_dir: Optional[str] = os.environ.get("MGTKDIR")
        if mgtk_dir is None:
            logger.warning(
                "MG-TK installation not found. \n"
                "MG-TK configuration will not be updated, as environmental "
                "variable MGTKDIR could not be found. \n"
                "Rerun get_gtdb.py configure on the system with MG-TK "
                "installed to complete configuration. \n"
                "Additionally, instructions on updating the configuration will "
                "be written to README.txt in the output "
                "directory and to the terminal at the end of the program."
            )
            return False

        logger.info("MGTK installation found in %s", mgtk_dir)
        x: str = prompt_input_set(
            "Do you want this script to automatically update MG-TK "
            "configuration to use the downloaded databases?",
            valid=["yes", "no", "cancel"],
            case=True
        )

        if x == "cancel":
            sys.exit()
        
        return x == "yes"

    def update_config(
        self,
        out_dir: Path
    ) -> None:
        """Attempt to automatically update MG-TK configuration to use new DB."""

        mgtk_dir_s: Optional[str] = str(self.get_mgtk_dir())
        if mgtk_dir_s is None:
            logger.error(
                "MG-TK installation not found. Please ensure environmental "
                "variable MGTKDIR contains path to MG-TK directory."
            )
            return
        mgtk_dir: Path = Path(mgtk_dir_s)

        logger.debug("Parse MG-TK config.txt at %s", mgtk_dir / "config.txt")
        mgtk_config: Dict[str, str] = self.parse_mgtk_config(
            mgtk_dir / "config.txt"
        )
        out_res: Path = out_dir.resolve()
        db_dir: Optional[str] = None
        use_abs: bool
        if "DBDir" in mgtk_config:
            db_dir = mgtk_config["DBDir"]
            # If databases are stored in the DBDir, use [DBDir] format, 
            # otherwise asbolute paths
            db_res: Path = Path(db_dir).resolve()
            logger.info("MG-TK DBDir: %s", db_dir)
            logger.debug("MG-TK DBDir resolved: %s", db_res)
            logger.debug("GTDB Database resolved: %s", out_res)
            if db_res in out_res.parents:
                use_abs = False
                logger.debug(
                    "Data in DBDir, using [DBDir] format."
                )
            else:
                logger.info(
                    "Downloaded database is not in DBDir, will use absolute "
                    "paths."
                )
                use_abs = True
        else:
            logger.info(
                "DBDir not found in MG-TK config.txt, will use absolute "
                "paths to database directories."
            )
            use_abs = True

        to_change: Dict[str, List[str]] = {
            'GTDBPath'      : ['markerGenes'],
            'GTDB_GTDB'     : [f'gtdb_r{self.version:.0f}_lineageGTDB.tab'],
            'GTDB_lnks'     : [f'gtdb_r{self.version:.0f}_clustering.tab'],
            'GTDBtk_DB'     : [f'gtdb'],
            'GTDBtk_mash'   : [f'gtdb', 'mashD'],
        }

        db_config: Path = mgtk_dir / "Mods" / "config_DBs.txt"
        i: int = 1
        while True:
            db_config_bup = db_config.with_suffix(f".bup{i}")
            if not db_config_bup.exists():
                break
            i += 1
        logger.info("Backing up current MG-TK config to %s", db_config_bup)
        shutil.copy(db_config, db_config_bup)

        mod_lines: List[str] = []
        with db_config.open("r") as f:
            for line in f:
                if len(line.strip()) < 1:
                    mod_lines.append(line)
                    continue
                if line.strip()[0] == '#':
                    mod_lines.append(line)
                    continue
                parts: List[str] = line.split("\t")
                if parts[0] in to_change:
                    # Update with new value
                    new_path: str = str(Path(out_res, *to_change[parts[0]]))
                    if not use_abs:
                        new_path = new_path.replace(str(db_res), "[DBDir]/")
                    mod_lines.append(
                        f"{parts[0]}\t{new_path}\t"
                        "#Updated by get_gtdb.py\n"
                    )
                    logger.debug("Updated %s to %s", parts[0], new_path)
                else:
                    mod_lines.append(line)

        logger.info("Writing updated config")
        with db_config.open("w") as f:
            f.writelines(mod_lines)

        # Attempt to install correct GTDB-tk
        ver_tup = self.required_gtdbtk(self.version)
        if ver_tup is None:
            logger.warning(
                "Cannot determine suitable GTDB-TK version, please check "
                "required version and install into environment MGTKgtdbtk "
                "manually."
            )
            return
        ver_min, ver_max = ver_tup
        logger.debug("Compatible GTDB-TK - Min: %s, Max %s", ver_min, ver_max)
        # Use the lowest compatible version if max version is current, as
        # no guarantee this script is up to date. Othewise, use higher
        # compatible version if fixed versions given.
        ver_use: str = (
            f"=={ver_min}" if ver_max == "Current" else f"=={ver_max}"
        )
        logger.info("Installing gtdbtk%s into environment MGTKgtdbtk", ver_use)
        cmd: str = (
            "micromamba install --name MGTKgtdbtk -c bioconda -c conda-forge "
            f"-y --channel-priority flexible gtdbtk{ver_use}"
        )
        logger.debug("Command: %s", cmd)
        subp: subprocess.CompletedProcess = subprocess.run(
            cmd, 
            shell=True,
            capture_output=True
        )
        if subp.returncode != 0:
            # Check if we got a n error due to fully received
            logger.error(
                "gtdbtk installation failed. Dumping output.")
            logger.error("Stderr")
            logger.error(subp.stderr)
            logger.error("Stdout")
            logger.error(subp.stdout)
            raise Exception(subp)

    def _processing_metadata(
            self,
            dest: Path,
            **kwargs
    ) -> None:
        """
        Output details of dates, versions and URLs to a JSON 'meta.json'
        
        :param dest: Output directory. Existing meta.json will be
        replaced without warning.
        """
        
        logger.info("Writing metdata to meta.json")

        meta: Dict = dict(
            version=self.version,
            urls=self.required_urls,
            date=str(datetime.now()),
            script_version=__version__,
            summary="""GTDB formatted for MG-TK using get_gtdb.py""",
        )

        with open(dest / "meta.json", "w") as f:
            json.dump(meta | kwargs, f, indent=4)

    def instructions(
        self,
        dest: Path,
        tk: bool
    ) -> str:
        """
        Write instructions for updating MG-TK and if needed getting GTDBtk data
        if it was not possible to do do automatically.

        :param dest: Extraction directory
        :param tk: Was GTDBtk data downloaded
        """

        instructions: str = (
            self._mgtk_instructions(dest)
            + "\n\n"
            + self._gtdbtk_instructions(dest, tk)
        )
        (dest / "README.txt").write_text(instructions)
        print(instructions)
        print("\nThese instructions were also written to:")
        print(str(dest / "README.txt"))
        return instructions

    def _mgtk_instructions(
            self,
            dest: Path
    ) -> str:
        """
        Compose instructions on update MG-TK config to use new database.
        
        String with the instructions will also be returned, so that it can be
        printed if desired.
        Currently this does not update config_DB.txt as a design decision,
        as don't want to have an unexpected and potentially impactful
        side effect if user just wanted to download and not replace, or not
        have in the directory the extracted to.

        :params dest: Destination directory.
        """


        mgtk_dir: Optional[str] = (
            "[Your MG-TK Dir]" if os.environ.get("MGTKDIR") is None else
            os.environ.get("MGTKDIR")
        )
        rec_ver: str = f"{self.version:.0f}"

        # Format a string with relevant parts
        instructions: str = (
            f"We recommend moving the output directory ({dest}) to your MGTK "
            "DBDir, and to a version specific subdirectory i.e. to\n\n"
            f"<DBDir>/MarkerG/GTDB_r{rec_ver}_MGTK\n\n"
            "To use this version in MG-TK, you must update config_DB.txt. This "
            f"is in {mgtk_dir}/Mods/config_DB.txt.\n"
            "Update the lines:\n\n"
            f"GTDBPath\t[DBDir]/MarkerG/GTDB_r{rec_ver}_MGTK/markerGenes/\n"
            f"GTDB_GTDB\t[DBDir]/MarkerG/GTDB_r{rec_ver}_MGTK/"
            f"gtdb_r{rec_ver}_lineageGTDB.tab\n"
            f"GTDB_lnks\t[DBDir]/MarkerG/GTDB_r{rec_ver}_MGTK/"
            f"gtdb_r{rec_ver}_clustering.tab\n\n"
            "Alternatively, you can run \n"
            f"get_gtdb.py configure -d {dest}\n"
            "on the system with MG-TK installed."
        )

        return instructions
    
    def _gtdbtk_instructions(
        self,
        dest: Path,
        tk: bool) -> str:
        """
        Instructions for GTDBtk data install or download.

        :param dest: Extraction directory
        :param tk: Was GTDBtk data downloaded
        """

        mgtk_dir: Optional[str] = (
            "[Your MG-TK Dir]" if os.environ.get("MGTKDIR") is None else
            os.environ.get("MGTKDIR")
        )
        rec_ver: str = f"{self.version:.0f}"
        ver: Optional[Tuple[str, str]] = self.required_gtdbtk(self.version)

        tk_inst: str = (
            (f"GTDBtk version between {ver[0]} and {ver[1]} is required for "
            "this release of GTDB. Please install this using micromamba into "
            "environment 'MGTKgtdbtk' i.e.\n\n"
            f"micromamba install --name MGTKgtdbtk gtbtk=={ver[0]}\n")
            if ver else
            ("We could not identify which version of GTDBtk is required "
            "for this release. Please check "
            "https://ecogenomics.github.io/GTDBTk/installing/"
            "index.html#gtdb-tk-reference-data and install the correct "
            "required version into the MGgtdbtk environment i.e.\n\n"
            "micromamba install --name MGTKgtdbtk gtdbtk==ver")
        )

        if tk:
            # Determine release specific subdirectory
            suffix: str
            try:
                subd: Path = next(x for x in dest.iterdir() if x.is_dir())
                suffix = f"/{subd.parts[-1]}"
            except StopIteration as e:
                suffix = ""
            tk_db: str = (
                "GTDBtk requires release specific databases. These have been "
                "download and extracted by this script. \n"
                "To use this version, you must update config_DBs.txt. This is "
                f"in {mgtk_dir}/Mods/config_DB.txt.\n"
                "Update the lines:\n\n"
                
                f"GTDBtk_DB\t[DBDir]/MarkerG/GTDB_r{rec_ver}_MGTK/"
                f"gtdb{suffix}\n"
                f"GTDBtk_mash\t[DBDir]/MarkerG/GTDB_r{rec_ver}_MGTK/gtdb"
                f"{suffix}/mashD"
            )
        else:
             tk_db: str = (
                "GTDBtk requires release specific databases. These were not "
                "downloaded by this script. \n"
                "Please see https://ecogenomics.github.io/GTDBTk/installing/"
                "index.html#gtdb-tk-reference-data for how to download this "
                "data. \n"
                "Additionally, GTDBtk has a script which will download the "
                "most current version, 'download-db.sh'. \n\n"
                "Once you have downloaded the correct version, you must update "
                "config_DB.txt. This is "
                f"in {mgtk_dir}/Mods/config_DB.txt.\n"
                "Update the lines:\n\n"
                
                f"GTDBtk_DB\t[DBDir]/MarkerG/GTDB_r{rec_ver}_MGTK/gtdbtk\n"
                f"GTDBtk_mash\t[DBDir]/MarkerG/GTDB_r{rec_ver}_MGTK/gtdbtk"
                "/mashD"
            )
        
        return tk_db + "\n\n" + tk_inst

    @staticmethod
    def get_mgtk_dir() -> Optional[Path]:
        """Get the MG-TK install path. Returns None if not found."""
        dir_s: Optional[str] = os.environ.get("MGTKDIR")
        if dir_s is None:
            return None
        return Path(dir_s)
    
    @staticmethod
    def parse_mgtk_config(cfg: Path) -> Dict[str, Optional[str]]:
        """Read MG-TK config files into a dictionary."""
        with cfg.open("r") as f:
            config: Dict[str, Optional[str]] = dict()
            for line in f:
                if len(line.strip()) < 1:
                    continue
                if line.strip()[0] == "#":
                    continue
                parts: List[str] = line.split("\t")
                k, v = parts[0], parts[1] if len(parts) > 1 else None
                config[k] = v
        return config

    @staticmethod
    def required_gtdbtk(version: float) -> Optional[Tuple[str, str]]:
        """
        Give GTDBtk requirements for this release of GTDB.

        Different GTDB release have different min and max versions of the
        GTDBtk software supported. This function will return which version
        is required, correct as of 12/02/2025. Where I have determined it,
        it will also return the URLs required to download the GTDBtk specific
        databases.

        The requirements are taken from 
        https://ecogenomics.github.io/GTDBTk/installing/index.html#gtdb-tk-reference-data
        one 12/02/2025.

        :param version: GTDB release
        :return: Tuple of min and max version
        """

        return (
            GTDBVersion.GTDBTK_VERSIONS[version]
            if version in GTDBVersion.GTDBTK_VERSIONS
            else None
        )

    def extract(
        self,
        dl_dir: Path = Path(os.getcwd()),
        dest_dir: Path = Path(os.getcwd()) / "output",
        tk: bool = False,
    ) -> None:
        """Format GTDB database for MG-TK.

        :param dl_dir: Path GTDB database was downloaded to
        :param dest_dir: Path to extract MG-TK formatted data to
        :param tk: Handle GTDBtk data
        """

        # Get a list of files which should have been downloaded
        files: Dict[str, Path] = self.download(
            download_dir=dl_dir,
            offline=True,
            download_tk=(False if tk == "skip" else tk)
        )

        # Extraction and formatting
        dest_dir.mkdir(parents=True, exist_ok=True)
        self._extract_markers(
            dest=dest_dir / "markerGenes",
            downloaded_files=files)

        self._format_taxonomy(dest=dest_dir, downloaded_files=files)
        if tk:
            self._extract_tk(
                dest=dest_dir / "gtdb",
                downloaded_files=files
            )

        self._processing_metadata(dest=dest_dir)

    @staticmethod
    def from_version(
        version: float,
        **kwargs
    ):
        """
        Get an object which will download a specific version of GTDB.

        This function checks the version and return a specific subclass if one
        exists to get a specific version. Otherwise, it creates a base class,
        which was developed to work on r220.

        :param version: Version of GTDB
        :param **kwargs: Passed to GTDBVersion.__init__
        """
        # Logic to return a version specific sublcass should go here when
        # I've identified any and implemented them.
        # Using default
        if version == 226:
            return GTDB226(version, **kwargs)
        if version != 220:
            logger.warning(
                "Using default download URLs and extraction methods\n"
                "This method was tested to work with r220, but are "
                f"untested for {version}."
            )
        return GTDBVersion(version, **kwargs)

class GTDB111(GTDBVersion):
    """This is an example of how you would override methods to implement
    downloads / extractions for new versions of GTDB in different formats.
    You should also add logic to select this class to the 
    GTDBVersion.from_version factory method."""

    def __init__(
            self
    ) -> None:
        """Example of how to set specific URLs for files."""
        super().__init__(
            version=111.0,
            required_urls=dict(
                bac_taxonomy="https://example.com/111/tax.tsv.gz",
                bac_markers="https://example.com/111/marker.tar.gz",
                bac_metadata="https://example.com/111/meta.tsv.gz",
                arc_taxonomy="https://example.com/111/tax.tsv.gz",
                arc_markers="https://example.com/111/marker.tar.gz",
                arc_metadata="https://example.com/111/meta.tsv.gz",
                tk_database="https://example.com/111/tkdb.tar.gz"
            )
        )
    
    def _extract_markers(self, dest, downloaded_files):
        """Implement version specific logic for extraction in here."""
        # You can access downloaded files by key, which returns local path
        raise NotImplementedError("Example class, no implementation")

    def _format_taxonomy(
            self, dest, downloaded_files, 
            taxonomy_headers = ["accession", "taxonomy"]
        ):
        """Implement version specific formatting of taxonomy here."""
        # If there were just different headers to look for in the source
        # taxonomy table, can just pass different arguments to the parent
        # method. Say 'accession' was called 'ncbi_accession' instead.
        super()._format_taxonomy(
            dest, downloaded_files,
            taxonomy_headers=["ncbi_accession", "taxonomy"])

class GTDB226(GTDBVersion):
    """Download and process v226 of GTDB for MG-TK
    
    The only significant change here is in how the taxonomic metadata is handled
    as marker genes can come from non-representative species."""

    def _make_gtdbmg_tax(
        self,
        dest: Path,
        tax: Path,
        lineage: Path
    ) -> None:
        """Make GTDBmg.tax allowing duplicates.
        
        MG-TK will attempt to automatically make a GTDBmg.tax file mapping 
        from accession to lineage. However, the method it uses will only keep 
        the last accession for the species. We will instead produce this file 
        retaining duplicates.

        :param dest: Directory to write GTDBmg.tax to
        :param tax: File with species to lineage mapping for GTDB
        :param lineage: File with lineage for each species
        """

        # Dictionary of species -> lineage lookup
        logger.info("Creating GTDBmg.tax with duplicates for species")
        logger.debug("Output to %s", dest / "GTDBmg.tax")
        lineage_map: Dict[str, Any] = dict()
        with lineage.open('rt') as f:
            for row in csv.reader(f, delimiter="\t"):
                species: str = row[0]
                lineage_s: str = ";".join(row[2:])
                lineage_map[species] = lineage_s

        with (
            tax.open('rt') as tax_in,
            (dest / "GTDBmg.tax").open('wt') as tax_out
        ):
            csv_in = csv.reader(tax_in, delimiter="\t")
            csv_out = csv.writer(tax_out, delimiter="\t")
            logger.debug("Write to %s", dest / "GTDBmg.tax")
            for row in csv_in:
                species: str = row[0]
                lineage_s = lineage_map[species]
                accession: str = row[1]
                csv_out.writerow(
                    [accession.strip(), lineage_s.strip()]
                )

    def _format_taxonomy(
            self,
            dest: Path,
            downloaded_files: Dict[str, Path],
            taxonomy_headers: List[str] = ["accession", "taxonomy"]
    ) -> None:
        """
        Make lineage and clustering tables in tab separated format. Not
        intended to be called directly, can be overriden for files which
        are not in the v220 format. This uses standard library functions to
        make it easier to use as part of generic install process.
        """

        logger.info("Producing taxonomy tables")
        logger.debug("Destination: %s", dest)

        logger.info("Concatenate source archaeal and bacterial tables")
        bac_tax, arc_tax = (
            downloaded_files['bac_taxonomy'], downloaded_files['arc_taxonomy'])
        bac_md, arc_md = (
            downloaded_files['bac_metadata'], downloaded_files['arc_metadata'])
        comb_md: Path = bac_md.parent / "combined_md.tsv.gz"
        comb_tax: Path = bac_md.parent / "combined_tax.tsv.gz"
        # Can't get process substitution working via subprocess, so having
        # to make some intermediate files (and don't want to handle in this 
        # script opening to append etc)
        cmd: str = (
            f"cat {bac_tax} {arc_tax} > {comb_tax} && "
            f"zcat {arc_md} | tail -n +2 | gzip > tmp_arc_md.tsv.gz && "
            f"cat {bac_md} tmp_arc_md.tsv.gz > {comb_md} && "
            "rm tmp_arc_md.tsv.gz"
        )
        logger.debug(f"Command: {cmd}")
        proc_res: subprocess.CompletedProcess = subprocess.run(
            cmd,
            capture_output=True,
            shell=True
        )
        proc_res.check_returncode()

        logger.info("Making species lineage")
        species: Set[str] = set()
        with (
            gzip.open(comb_tax, 'rt') as tbl,
            open(
                dest / f"gtdb_r{self.version:.0f}_lineageGTDB.tab",
                "w") as lin_o
        ):
            lineage_write = csv.writer(
                lin_o, delimiter="\t", quoting=csv.QUOTE_NONE)
            for row in csv.DictReader(
                tbl, delimiter="\t", fieldnames=taxonomy_headers):
                accession, lineage = row['accession'], row['taxonomy']
                lineage_parts: List[str] = [
                    x.strip() for x in lineage.split(";")]
                # Only add a species once
                if lineage_parts[-1] in species:
                    continue
                species.add(lineage_parts[-1])
                lineage_write.writerow(
                    [lineage_parts[-1].strip(), 0, *lineage_parts]
                )

        logger.info("Making species to accession mapping")
        logger.info(
            "For v226, we output species for all accessions as some genes "
            "are not from representatives"
        )
        with (
            gzip.open(comb_md, 'rt') as tbl,
            open(
                dest / f"gtdb_r{self.version:.0f}_clustering.tab",
                "w") as clust_o
        ):
            lineage_write = csv.writer(
                clust_o, delimiter="\t", quoting=csv.QUOTE_NONE)
            for row in csv.DictReader(tbl, delimiter="\t"):
                accession: str = row['accession']
                # if accession.strip() in accessions:
                if True:
                    tax: str = row['gtdb_taxonomy']
                    species_s = tax.split(";")[-1]
                    lineage_write.writerow(
                        [species_s.strip(), accession.strip()])
        
        self._make_gtdbmg_tax(
            dest / "markerGenes",
            lineage=dest / f"gtdb_r{self.version:.0f}_lineageGTDB.tab",
            tax=dest / f"gtdb_r{self.version:.0f}_clustering.tab"
        )

        comb_md.unlink()
        comb_tax.unlink()

        logger.info("Taxonomy tables finished")


def is_url_dir(src: str) -> bool:
    return src[-1] == "/"

def mkdir_parents(dest: str) -> None:
    dest_pth: Path = Path(dest)
    logger.debug("Making directory %s", dest_pth)
    if "." not in dest_pth.name:
        dest_pth.mkdir(parents=True, exist_ok=True)
    else:
        dest_pth.parents[0].mkdir(parents=True, exist_ok=True)

def cmd_wget(src: str, dest: str) -> str:
    """Compose a command to fetch a file or directory using wget."""
    mkdir_parents(dest)
    if is_url_dir(src):
        # Number of directories to trim from structure
        return (
            f"wget -r -np -nH -nd -R 'index.html*' --level=1 -nc "
            f"-c {src} -P {dest}"
        )
    else:
        return f"wget -c {src} -O {dest}"

def cmd_curl(src: str, dest: str) -> str:
    """Compose a command to fetch a file or directory using curl."""
    mkdir_parents(dest)
    if is_url_dir(src):
        raise NotImplementedError(
            "Directory download is not supported by curl. "
            "Please install wget, or download full GTDBtk package."
        )
    else:
        return f"curl -o {dest} -C - {src}"

def cmd_test(src: str, dest: str) -> str:
    """Compose a command to create a dummy file for a download.
    This uses some preconstructed dummy files which are in 
    /helpers/install, and if using test mode expects to be run in
    this helpers/install."""
    mkdir_parents(dest)
    dummy_root: Path = Path("get_gtdb/")
    if "marker_genes" in src:
        if "ar" in src:
            return f"cp {dummy_root / 'arc_markers.tar.gz'} {dest}"
        elif "bac" in src:
            return f"cp {dummy_root / 'bac_markers.tar.gz'} {dest}"
    elif "taxonomy"in src:
        if "ar" in src:
            return f"cp {dummy_root / 'arc_taxonomy.tsv.gz'} {dest}"
        elif "bac" in src:
            return f"cp {dummy_root / 'bac_taxonomy.tsv.gz'} {dest}"
    elif "metadata" in src:
        if "ar" in src:
            return f"cp {dummy_root / 'arc_metadata.tsv.gz'} {dest}"
        elif "bac" in src:
            return f"cp {dummy_root / 'bac_metadata.tsv.gz'} {dest}"
    elif "split_package" in src:
        # GTDB split package
        # Make 3 fake split parts

        parts = [
            (
                f'cp {dummy_root}/gtdbtk_dummy.tar.gz.part_{p} '
                f'{dest}/gtdbtk_dummy.tar.gz.part_{p}'
            )
            for p in 
            ['aa', 'ab', 'ac']
        ]
        return " && ".join(parts)
    elif "full_package" in src:
        return f"touch {dest}"
    logger.warning("Defaulting to empty file for %s in %s", src, dest)
    return f"touch {dest}"

def prompt_input_set(prompt: str, valid: List[str], case: bool = False) -> str:
    """Prompt user for input from a set of options. If not case sensitive
    returned value is always lower case. Will try partial matching stem of
    word."""
    valid = [i.lower() if case else i for i in valid]
    prompt_form: str = f"{prompt} ({'/'.join(valid)}): "
    while True: 
        x: str = input(prompt_form)
        x = x.lower() if case else x
        x_len: int = len(x)
        matches: List[str] = [i for i in valid if i[:x_len] == x]
        if len(matches) == 1:
            logger.debug("User entered %s, selected %s", x, matches[0])
            return matches[0]
        logger.warning("%s not valid, must be one of %s", x, str(valid))

def greet():
    """Print a nicely formatted welcome message"""
    print("""
\033[96mget_gtdb\033[0m
\033[96m────────\033[0m""")

def add_args(
    parser,
    *args
):
    """Add an argument to a command parser. arg should be a tuple of positional
    then keyword args."""
    for arg in args:
        parser.add_argument(
            *arg[0],
            **arg[1]
        )

def configure_logging(debug: bool = False):
    """Configure the logging level based on command line flag `debug`"""
    if debug:
        logger.setLevel(logging.DEBUG)

def meta_from_dir(dir: Path) -> Dict[str, Any]:
    """Get settings from a directory."""
    json_pth: Path = dir / "meta.json"
    if not json_pth.exists():
        logger.error(
            "No meta.json found in %s, unable to determine version. Please "
            "provide a directory which was created by this tool."
        )
        raise FileNotFoundError()
    with json_pth.open("rt") as f:
        meta_d = json.load(f)
    return meta_d

def cli_download(args):
    configure_logging(args.debug)
    greet()
    downloader: GTDBVersion = GTDBVersion.from_version(
        args.version,
        test=args.test
    )
    downloader.download(
        download_dir = args.tmp,
        download_tk = args.tk
    )

def cli_extract(args):
    configure_logging(args.debug)
    greet()
    meta = meta_from_dir(args.tmp)
    version, tk = meta['version'], meta['tk']
    downloader: GTDBVersion = GTDBVersion.from_version(version)
    downloader.extract(
        dl_dir = args.tmp,
        dest_dir = args.dest,
        tk = tk
    )

def cli_configure(args):
    configure_logging(args.debug)
    greet()
    meta = meta_from_dir(args.dest)
    version = meta['version']
    downloader: GTDBVersion = GTDBVersion.from_version(version)
    downloader.update_config(
        out_dir = args.dest
    )

def cli_all(args):
    """Run all steps - download, extract, configure."""
    configure_logging(args.debug)
    greet()
    version = args.version
    downloader: GTDBVersion = GTDBVersion.from_version(
        version,
        test=args.test
    )
    # Ask if want to update config
    config: bool = GTDBVersion.ask_config_update()
    downloader.download(
        download_dir = args.tmp,
        download_tk = args.tk
    )
    downloader.extract(
        dl_dir = args.tmp,
        dest_dir = args.dest,
        tk = args.tk
    )
    if config:
        downloader.update_config(
            out_dir = args.dest
        )
    else:
        logger.warning("Did not update MG-TK configuration update.")
        logger.warning(
            "To update MG-TK to use new database version, either run "
            "'get_gtdb.py configure -d %s' on the system with MG-TK installed "
            "or follow the instructions in %s.",
            str(args.dest),
            str(args.dest / "README.txt")
        )
        downloader.instructions(
            dest=args.dest,
            tk=args.tk
        )

def main():
    """Command line interface"""


    parser = argparse.ArgumentParser(
        prog="get_gtdb.py",
        description="""
            Download and format versions of the GTDB database for use in 
            the pipeline MG-TK. This program expects GNU tar and either wget
            or curl to be available.
        """,
        epilog="""
            This is separated into three steps, which can be run separately as
            subcommands 'download', 'extract', and 'configure'. They can also
            be all run at once using subcommand 'all'.

            This program expects GNU tar and either wget or curl to be 
            available. Downloads try to resume if they failed, and should be 
            skipped for existing files (of the same size) when the script is
            rerun.
            
            If the system you want to put the created database on lacks
            internet access, use 'download' subcommand to download files, then
            transfer to target system, then use 'extract' to process files.
            
            'configure' will update MG-TK configuration to use an extracted
            database, and install an appropriate version of gtdbtk. If your
            group is sharing copies of the database, once it has been extracted,
            others should be able to configure MG-TK to use it using the.
            'configure' sucommand.
        """
    )

    # Define some arguments which will be used across multiple subparsers
    arg_dest = (
        ["-d", "--dest"],
        dict(
            type=Path,
            required=False,
            default=Path(os.getcwd()) / "output",
            help=(
                "Directory for MG-TK formatted database. When extracting, "
                "contents will be overwritten."
            )
        )
    )
    arg_version = (
        ["-v", "--version"],
        dict(
            type=float,
            required=True,
            help="Version of GTDB to download. Can include minor versions (202.1)"
        )
    )
    arg_temp = (
        ["-t", "--tmp"],
        dict(
            type=Path,
            required=False,
            default=Path(os.getcwd()),
            help=(
                "Temporary directory to download files from GTDB to."
            )
        )
    )
    arg_tk = (
        ["--tk"],
        dict(
            dest="tk",
            type=str,
            choices=["split", "full", "skip"],
            default="split",
            help=(
                "Method used to download the GTDBtk database. This is very large "
                "(~110GB for r220). By default this will be downloaded in parts "
                "which are then concatenated ('split'). You can instead download "
                "the full file ('full'). This can also be be skipped using 'skip' "
                "if you already have this data or will acquire it another way."
            )
        )
    )
    arg_test = (
        ["--test"],
        dict(
            dest="test",
            action="store_true",
            help=(
                "Do not download any files, instead create dummy files in the"
                "destination. Included to allow test of script logic without "
                "needing to download and extract large files."
            )
        )
    )



    # Always used arguments
    parser.add_argument(
        "--debug",
        action="store_true",
        help=("Produce debug messages. This will output any shell commands "
              "being run.")
    )

    # Subcommand parsers
    subparsers = parser.add_subparsers(help="Subcommand help")

    # Download subcommand
    parser_download = subparsers.add_parser(
        "download",
        help="Download GTDB and GTDBtk databases.",
        description="Download GTDB and GTDBtk databases."
    )
    parser_download.set_defaults(func=cli_download)
    add_args(parser_download, arg_temp, arg_version, arg_tk, arg_test)

    # Extract subcommand
    parser_extract = subparsers.add_parser(
        "extract",
        help=(
            "Extract downloaded GTDB and GTDBtk databases to the format"
            "and structure required by MG-TK."
        ),
        description=(
            "Extract downloaded GTDB and GTDBtk databases to the format"
            "and structure required by MG-TK."
        )
    )
    parser_extract.set_defaults(func=cli_extract)
    add_args(parser_extract, arg_temp, arg_dest)

    # Configure subcommand
    parser_configure = subparsers.add_parser(
        "configure",
        help=(
            "Update MG-TK configuration and install correct version of GTDBtk "
            "to use a downloaded and extracted database version."
        ),
        description=(
            "Update MG-TK configuration and install correct version of GTDBtk "
            "to use a downloaded and extracted database version. Will require "
            "internet access for MG-TK download. Should be run by the user "
            "with MG-TK installed."
        )
    )
    parser_configure.set_defaults(func=cli_configure)
    add_args(parser_configure, arg_dest)

    # All subcommand
    parser_all = subparsers.add_parser(
        "all",
        help=(
            "Download and extract GTDB and GTDBtk, and configure MG-TK to "
            "use download version. Equivalent to runining download, extract, "
            "then configure subcommands. Will required internet access for "
            "download, and GTDBtk install."
        )
    )
    parser_all.set_defaults(func=cli_all)
    add_args(parser_all, arg_temp, arg_dest, arg_version, arg_tk, arg_test)

    args = parser.parse_args()
    if "func" in args:
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
