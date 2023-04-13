from pyrosetta import *
from pyrosetta.toolbox import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta.protocols.simple_moves import *
from pyrosetta.rosetta.protocols.minimization_packing import *
from Bio import AlignIO
import numpy as np
import os,sys,re,random,shutil,subprocess
init()

class CALCSCORE():
	def __init__(self,pose,outputfile):
		self.scorefxn = ScoreFunction()
		self.scorefxn = get_fa_scorefxn()

		self.refine = Pose()
		self.refine.assign(pose)

		#Definition of energy minimize protocol
		#Start
		self.kT = 1.0
		self.cycles = 100
		self.movemap = MoveMap()
		self.movemap.set_bb(False)
		#movemap.set_bb(True)

		self.minmover = MinMover()
		self.minmover.movemap(self.movemap)
		self.minmover.score_function(self.scorefxn)
		
		self.to_pack = standard_packer_task(pose)
		self.to_pack.restrict_to_repacking()    # prevents design, packing only
		self.to_pack.or_include_current(True)    # considers the original sidechains
		self.packmover = PackRotamersMover(self.scorefxn, self.to_pack)

		self.combined_mover = SequenceMover()
		self.combined_mover.add_mover(self.minmover)
		self.combined_mover.add_mover(self.packmover)

		self.mc = MonteCarlo(pose,self.scorefxn,self.kT)
		self.trial = TrialMover(self.combined_mover,self.mc)
		self.refinement = RepeatMover(self.trial,self.cycles)
		#End of EM protocol.

		self.refine = Pose()
		self.refine.assign(pose)

		self.refinement.apply(self.refine)

		self.score_refine = self.scorefxn(self.refine)
		self.score_init = self.scorefxn(pose)

		self.refine.dump_pdb("%(outputfile)s.pdb"%vars())

		print(f"The rosetta score prior to the energy minimization: {self.score_init}\n")
		print(f"The rosetta score after the energy minimization: {self.score_refine}\n")

while len(sys.argv)>1:
	option = sys.argv[1]
	del sys.argv[1]
	if option == "-PDB":
		inputpdb = sys.argv[1]
		del sys.argv[1]
	elif option == '-OUTPUT':
		outputfile = sys.argv[1]
		del sys.argv[1]

pose = Pose()
pose = pose_from_pdb(inputpdb)
CALCSCORE(pose,outputfile)





