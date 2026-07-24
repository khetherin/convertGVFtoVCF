import os
import yaml


class ProjectPaths:
    """ The responsibility of this class is to manage paths from a config file.
    It resolves the full path for the files in the config.
    """
    def __init__(self, config_path="config.yaml"):
        cfg_path = config_path if config_path is not None else "config.yaml"
        self.root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) # TOP LEVEL
        self.package_dir = os.path.dirname(os.path.abspath(__file__)) # project directory of convertGVFtoVCF
        self.full_config_path = os.path.join(self.package_dir, "etc", cfg_path)
        # reading in
        try:
            with open(self.full_config_path, 'r') as f:
                # data = yaml.safe_load(f) or {}
                data = yaml.safe_load(os.path.expandvars(f.read())) or {}
        except FileNotFoundError:
            raise FileNotFoundError(f"Required config missing at: {self.full_config_path}")
        # internal paths section of config
        paths = data.get('paths', {})
        self.etc_dir = os.path.join(self.package_dir, paths.get("etc_folder", ""))
        self.eva_schema_path = os.path.join(self.package_dir, paths.get("eva_schema_file", ""))
        self.dgva_schema_path = os.path.join(self.package_dir, paths.get("dgva_schema_file", ""))
        self.test_dir = os.path.normpath(os.path.join(self.package_dir, paths.get("test_folder", "")))
        # external paths section of config
        self.assembly_paths = data.get('assembly_paths', {})
        self.assembly_report_paths = data.get('assembly_report_paths', {})

    def get_assembly_path(self, assembly_input):
        """ Gets the assembly path
        :params assembly_input: assembly name in config.yaml e.g. GRCh38 or path
        :return full_path: assembly full path
        """
        if not assembly_input:
            return None
        if assembly_input in self.assembly_paths:
            # found the assembly in config.yaml
            target_path = self.assembly_paths[assembly_input]
        else:
            target_path = assembly_input
        if "$" in target_path:
            # to resolve ${REF_PATH} in config.yaml
            target_path = os.path.expandvars(target_path)
        # get full path
        if os.path.isabs(target_path):
            full_path = target_path
        else:
            full_path = os.path.join(self.root_dir, target_path)

        return os.path.normpath(full_path)