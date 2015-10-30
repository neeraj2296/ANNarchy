"""

    ProfileGenerator.py

    This file is part of ANNarchy.

    Copyright (C) 2013-2016  Julien Vitay <julien.vitay@gmail.com>,
    Helge Uelo Dinkelbach <helge.dinkelbach@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ANNarchy is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
import ANNarchy.core.Global as Global
from .Template.ProfileTemplate import profile_template

class ProfileGenerator(object):
    """
    Extent the generated code by profiling annotations.
    """
    def __init__(self, annarchy_dir):
        """
        Initialize ProfileGenerator.
        """
        self.annarchy_dir = annarchy_dir

    def generate(self):
        """
        Generate Profiling class code, called from Generator instance.
        """
        # Generate header for profiling
        with open(self.annarchy_dir+'/generate/Profiling.h', 'w') as ofile:
            ofile.write(self._generate_header())

    def generate_init_population(self, pop):
        """
        Generate initialization code for population
        """
        declare = """
    Measurement* measure_step;
"""
        init = """        // Profiling
        measure_step = Profiling::get_instance()->register_function("pop", "%(name)s", "step");
""" % { 'name': pop.name }

        return declare, init

    def generate_init_projection(self, proj):
        """
        Generate initialization code for projection
        """
        declare = """
    Measurement* measure_psp;
    Measurement* measure_step;
"""
        init = """        // Profiling
        measure_psp = Profiling::get_instance()->register_function("proj", "proj%(id_proj)s", "psp");
        measure_step = Profiling::get_instance()->register_function("proj", "proj%(id_proj)s", "step");
""" % { 'id_proj': proj.id }

        return declare, init

    def annotate_computesum_rate_omp(self, proj, code):
        """
        annotate the computesum compuation code
        """
        prof_begin = profile_template['compute_psp']['before']
        prof_end = profile_template['compute_psp']['after']

        prof_code = """
        // first run, measuring average time
        %(prof_begin)s
%(code)s
        %(prof_end)s
""" % {'code': code,
       'prof_begin': prof_begin,
       'prof_end': prof_end
       }
        return prof_code

    def annotate_computesum_spiking_omp(self, proj, code):
        """
        annotate the computesum compuation code
        """
        prof_begin = profile_template['compute_psp']['before'] % { 'name': 'proj'+str(proj.id) }
        prof_end = profile_template['compute_psp']['after'] % { 'name': 'proj'+str(proj.id) }

        prof_code = """
        // first run, measuring average time
        %(prof_begin)s
%(code)s
        %(prof_end)s
""" % {'code': code,
       'prof_begin': prof_begin,
       'prof_end': prof_end
       }
        return prof_code

    def annotate_update_synapse_omp(self, proj, code):
        """
        annotate the update synapse code, generated by ProjectionGenerator.update_synapse()
        """
        prof_begin = profile_template['update_synapse']['before']
        prof_end = profile_template['update_synapse']['after']

        prof_code = """
// first run, measuring average time
%(prof_begin)s
%(code)s
%(prof_end)s
""" % {'code': code,
       'prof_begin': prof_begin,
       'prof_end': prof_end
       }

        return prof_code

    def annotate_update_neuron_omp(self, pop, code):
        """
        annotate the update neuron code
        """
        prof_begin = profile_template['update_neuron']['before'] % { 'name': pop.name }
        prof_end = profile_template['update_neuron']['after'] % { 'name': pop.name }

        prof_code = """
        // first run, measuring average time
        %(prof_begin)s
%(code)s
        %(prof_end)s
""" % {'code': code,
       'prof_begin': prof_begin,
       'prof_end': prof_end
       }
        return prof_code

    def _generate_header(self):
        """
        generate Profiling.h
        """
        from .Template.ProfileTemplate import profile_header

        if Global.config["paradigm"] == "openmp":
            config_xml = """
        _out_file << "  <config>" << std::endl;
        _out_file << "    <paradigm>%(paradigm)s</paradigm>" << std::endl;
        _out_file << "    <num_threads>%(num_threads)s</num_threads>" << std::endl;
        _out_file << "  </config>" << std::endl;
        """ % { 'paradigm': Global.config["paradigm"], 'num_threads': Global.config["num_threads"]}
            config = Global.config["paradigm"] + '_'  + str(Global.config["num_threads"]) + 'threads'
            return profile_header % { 'config': config, 'config_xml': config_xml }
        else:
            return ""
