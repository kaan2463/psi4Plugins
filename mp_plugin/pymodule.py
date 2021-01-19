#
# @BEGIN LICENSE
#
# mp_plugin by Psi4 Developer, a plugin to:
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

import psi4
from psi4 import core
import psi4.driver.p4util as p4util
from psi4.driver.procrouting import proc_util

def run_mp_plugin(name, **kwargs):
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = psi4.driver.scf_helper(name, use_c1=True, **kwargs)  # C1 certified
    else:
        if ref_wfn.molecule().schoenflies_symbol() != 'c1':
            raise ValidationError("""  CBASIC2 does not make use of molecular symmetry: """
                                  """reference wavefunction must be C1.\n""")
    mp_plugin_wfn = psi4.core.plugin('mp_plugin.so', ref_wfn)

    return mp_plugin_wfn


psi4.driver.procedures['energy']['mp_plugin'] = run_mp_plugin
