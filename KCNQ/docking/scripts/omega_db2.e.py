#!/usr/bin/env python2.7

from __future__ import division


VERSION = "0.1.0"
DESCRIPTION = """Generate an ensemble of low energy conformations for a substance using Open Eye Omega as input for ddb2 format"""
USAGE = """

Usage:
    python omega_db2.e.py [options] <substance>.mol2

Output:
    ring_count
      a file containing a single line with the number of rings in the substance

    <substance>.<ring_index>.db2in.mol2
      where <ring_index> 0 if there are no rings or in [1-<ring_count>] otherwise
      containing a set of conformations for the input substance.

    file exit code is the number of substances omega failed to load

Details:
    The size and diversity of the set of conformations is a major contribution to the
    depth of sampling in Dock. Sampling more deeply allows findinging poses that fit
    better into binding sites--leading to lower enthalpy. However, on average, sampling
    deeper finds conformations that are in narrower energy wells--leading to lower entropy.
    Additionally, sampling deeper leads to increased computational cost per-compound.
    Since the number of conformations increases exponentially with the number of rotatable bonds,
    some care must be taken in controling how compounds are generated.

    Options for OEOmega are described here:
       https://docs.eyesopen.com/toolkits/python/omegatk/OEConfGenClasses/OEOmegaOptions.html
       https://docs.eyesopen.com/applications/omega/omega/omega_opt_params.html

    To balance configurability and convention Omega parameters are set
    by the following methods with descending priority:

       positional arguments
          numHs               (second positional argument after the input .mol2 file)
       Command line arguments
          fixrms
          rms
          maxconfs
          ewindow
       Environment variables
          OMEGA_ENERGY_WINDOW (alias for ewiondow)
          OMEGA_MAX_CONFS     (alias for maxconfs)
       hard coded non-standard options (see code below)
       omega default values

   If you find there is an option you would like to set for your application, the recommendation is to
   add it to the command line arguments and pass in the non-default value as a flag.

Examples:
  To increase sampling
    python omega_db2.e.py \\
       --fixrms 0.1 \\
       --rms 0.01 \\
       --maxconfs 1000 \\
       --zero_ring_maxconfs 1000 \\
       <substance>.mol2
"""
import os
import sys
import logging
import optparse
from openeye.oechem import *
from openeye.oeomega import *


def set_defaults(omega):
    # File Parameters
    omega.SetCommentEnergy(True)
    #omega.SetIncludeInput(False)
    omega.SetIncludeInput(True)
    omega.SetRotorOffset(False) #Enkhee
    omega.SetSDEnergy(False)
    omega.SetWarts(True)
    # 3D Construction Parameters
    omega.SetBuildForceField('mmff94s')
    omega.SetCanonOrder(True)
    omega.SetFixDeleteH(True)
    omega.SetDielectric(1.0)
    omega.SetExponent(1.0)
    omega.SetFixRMS(0.15)
    #omega.SetFixRMS(0.01)# jklyu modify this to keep more conformors
    omega.SetFromCT(False)
    omega.SetFixMaxMatch(1)
    omega.SetFixUniqueMatch(True)
    # Structure Enumeration Parameters
    omega.SetEnumNitrogen(False)
    omega.SetEnumRing(False)#Enkhee  # TS & MK 20160524 (from False) (Improves ring pucker sampling)
    # Torsion Driving Parameters
    omega.SetEnergyWindow(200)  # JJI 20160329  # TS & MK 20160524 (from 6)
    omega.SetMaxConfs(12)
    omega.SetMaxRotors(-1)
    omega.SetMaxSearchTime(120.0)
    omega.SetRangeIncrement(5)
    omega.SetRMSThreshold(0.50)  # JJI 20160329 # TS & MK 20160524 (from .5)
    #omega.SetRMSThreshold(0.01)  # jklyu modify this to keep more conformors
    omega.SetSearchForceField('mmff94s')
    
    omega.SetTorsionDrive(True)
    #omega.SetTorsionDrive(False)
    # Stereochemsitry
    #omega.SetStrictStereo(False)
    omega.SetStrictAtomTypes(False)

def write_molecule(outfile, mol):
    ofs = oemolostream(outfile)
    OEWriteMolecule(ofs, mol)
    ofs.close()

def generate_conformations(omega, infile, zero_ring_maxconfs):

    inroot, inext = os.path.splitext(infile)
    
    ## read read in file
    mol = OEMol()
    ifs = oemolistream(infile)
    OEReadMolecule(ifs, mol)
    ifs.close()

    #### Write out file reporting the nubmer of rings in the molecule
    OEDetermineComponents(mol)
    count, ringlist = OEDetermineRingSystems(mol)
    rcf = open('ring_count', 'w')
    rcf.write('%d\n' % count)
    rcf.close()

    ### Write out conformations for each ring system
    outfiles = []
    fail_count = 0
    print ('energy: ', omega.GetEnergyWindow())
    print ('conf: ', omega.GetMaxConfs() )
    if count == 0:
        omega.SetMaxConfs(zero_ring_maxconfs)
        outfile = "%s.%d.db2in.mol2" % (inroot, count) 
        outfiles.append(outfile)
        if omega(mol):
            write_molecule(outfile, mol)
        else:
            fail_count +=1 
    else:
        pred = OEPartPredAtom(ringlist)
        for i in xrange(1, count+1):
            pred.SelectPart(i)
            outfile = "%s.%d.db2in.mol2" % (inroot, i) 
            outfiles.append(outfile)
            molcopy = OEMol(mol)
            if omega(molcopy, pred):
                write_molecule(outfile, molcopy)
            else:
                fail_count +=1
                
    print outfiles

    return fail_count


def main(argv):
    parser = optparse.OptionParser(version=VERSION, description=DESCRIPTION, usage=USAGE)
    parser.add_option(
        "--fixrms", type="float", action="store", dest="fixrms", default=.15,
        help="""Fixfile fragments taken from crystallographic sources may differ in
                their geometry relative to optimal MMFF94
                geometries. OMEGA attempts to superimpose built
                structures onto fixfile fragments and, if the geometry
                differs too greatly, OMEGA considers the superposition
                a poor match and will fail to build a structure using
                the fixfile. This flag can be used to loosen the
                default RMS superposition criteria to allow suboptimal
                superpositions to succeed in spite of the poor
                geometric complementarity. [default = 0.15]""")
    parser.add_option(
        "--rms", type="float", action="store", dest="rms", default=0.5,
        help="Duplication threshold. Default value is 0.5 A")
    parser.add_option(
        "--maxconfs", type="int", action="store", dest="maxconfs", default=200,
        help="""sets the maximum number of conformations to be
                generated. Conformers are assembled in energy sorted
                order. As a special case, setting -maxconfs 0 will
                result in OMEGA skipping the duplicate removal step
                and it will write all generated conformers to the
                output file. Note that this implies -rms 0 is also
                used. [default = 200]""")
    parser.add_option(
        "--zero_ring_maxconfs", type="int", action="store", dest="zero_ring_maxconfs", default=30,
        help="""Set the maximum number of conformations for compounds
                that have zero ring systems""")
    parser.add_option(
        "--ewindow", type="float", action="store", dest="ewindow", default=12,
        help="""sets the energy window in kcal/mol used as an accept or reject
              criteria for conformers. Any conformer that has a
              calculated strain energy less than the sum of the energy
              window and the energy of the global minimum conformer
              will be accepted. Conformers with strain energies above
              this threshold are rejected. [default=12.0]""")
    parser.add_option(
        "--numHs", type="int", action="store", dest="numHs", default="-1",
        help="""Report number of rotatable hydrogens and use to scale max number of conformers
                if 2-3 --> reduce maxconfs by a factor of 3
                if 4-5 --> reduce maxconfs by a factor of 30
                if > 5 refuse to build conformations with > 5 rotatable bonds""")

    options, args = parser.parse_args()

    ### initialize omega
    omega = OEOmega()
    set_defaults(omega)

    ## Set the energy window
    if os.environ.get('OMEGA_ENERGY_WINDOW', '').strip() != '':
        ewindow = int(os.environ['OMEGA_ENERGY_WINDOW'])
        logging.warn("OMEGA_ENERGY_WINDOW from environment variable is set to {}.".format(ewindow, ewindow))
    else:
        ewindow = options.ewindow
    omega.SetEnergyWindow(ewindow)
        
    
    ## set the maximum number of conformations
    if os.environ.get('OMEGA_MAX_CONFS', '').strip() != '':
        maxconfs = int(os.environ['OMEGA_MAX_CONFS'])
        logging.warn("OMEGA_MAX_CONFS from environment variable is set to {}.".format(maxconfs, maxconfs))
    else:
        maxconfs = options.maxconfs

    # optionally scale the maximum number of conformations depending on the number of rotatable hydrogens
    if len(args) > 1:
        numHs = args[1]
        logging.warn(
            "Attempting to set numHs via positional argument, Please use --numHs={} instead.".format(numHs))
        assert numHs.isdigit()
        numHs=int(numHs)
    else:
        numHs = options.numHs
    if numHs > 5:
        print("Refusing to build conformations with > 5 rotatable hydrogens")
        sys.exit(-1)
    if numHs >= 4:
        logging.warn("4-5 Rotatable hydrogens  reported. Reducing confs by a factor of 30")
        maxconfs = maxconfs // 30
    elif numHs >= 2:
        logging.warn("2-3 Rotatable hydrogens  reported. Reducing confs by a factor of 3")
        maxconfs = numconfs // 3
    omega.SetMaxConfs(maxconfs)
    
    
    ### Initialize input file
    infile = args[0]

    fail_count = generate_conformations(omega, infile, zero_ring_maxconfs=options.zero_ring_maxconfs)
    
    sys.exit(fail_count)

if __name__ == "__main__":
    main(sys.argv)
