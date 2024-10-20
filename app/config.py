import os
import configparser

class Config:
    def __init__(self, path='./'):
        self._config_dir = path
        self._config_name = 'quickspec.ini' 
        self._config_path = os.path.join(self._config_dir, self._config_name)

        self.config = configparser.ConfigParser(inline_comment_prefixes=('#', ';')) 

        if not os.path.isfile(self._config_path):
            self._new_config()

        self.config.read(self._config_path)

    def save(self):
        with open(self._config_path, 'w') as configfile:
            self.config.write(configfile)

    def _check_section(self, section):
        if not self.config.has_section(section):
            self.config.add_section(section)

    def _new_config(self):
        self._check_section('window')
        self.config['window']['geometry'] = '600x672+330+111'
        self.save()
