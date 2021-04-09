#!/usr/bin/env python3
from lib import cc as virus
from vaccine.load import dat as vaccine
from corona import corona
from viralDNA import viralDNA

print("\n\nWelcome to: Reverse engineering the coronavirus\n---------------------------")

# Allows a maximum of 25 gene snippit lengths to be provided
# These are later passed to the virualDNA prototype
vaccine_i = 21562

regions = {
	"untranslated_region_f" : 265, 
	"orf1a_i" : 266 ,
	"orf1a_f" : 13483, 
	"orf1b_i" : 13468, 
	"orf1b_f" : 21555, 
	"spike_gp_i" : 21563, 
	"spike_gp_f" : 25384, 
	"orf3a_i" : 25393, 
	"orf3a_f" : 26220, 
	"envelope_p_i" : 26245, 
	"envelope_p_f" : 26472, 
	"membrane_gp_i" : 26523, 
	"membrane_gp_f" : 27191, 
	"orf6_i" : 27202,
	"orf6_f" : 27387,
	"orf7a_i" : 27394,
	"orf7a_f" : 27759,
	"orf7b_i" : 27756,
	"orf7b_f": 27887, 
	"orf8_i" : 27894, 
	"orf8_f" : 28259, 
	"n_p_i" : 28274, 
	"n_p_f" : 29533, 
	"orf10_i" : 29558, 
	"orf10_f" : 29674, 
	"new_region_i" : 0, 
	"new_region_f" : 0}

# User has 10 chances to change the regions 
for i in range(10):
	print("1 - Add a new gene region\n2 - Change a gene region value\n3 - Add the vaccine start index\n4 - Translate/GO") 

	# User wants to add a new gene region
	user_action = int(input())
	if (user_action == 1): 
		print("Press 0 to stop. Enter the start index of the new region: ?")
		new_region_i = int(input())

		if (new_region_i != 0):
			print("Press 0 to stop. Enter the final index of the new region: ?")
			new_region_f = int(input())

	elif (user_action == 2):
		# for region in region_names:
		print("Which region do you want to change?\n")
		print("Type 0 to go back to Main Menu. Type -1 to skip. Choose a menu option:")
		print("1 - Start position of  orf1a (Default = 266)\n")
		print("2 - End position of orf1a (Default = 13483)\n")
		print("3 - Start position of orf1b (Default = 13468)\n")
		print("4 - End position of orf1b (Default = 21555)\n")
		print("5 - Start position of spike glycoprotein (Default = 21563)\n")
		print("6 - End position of spike glycoprotein (Default = 25384)\n")
		print("7 - Start position of orf3a (Default = 25393)\n")
		print("8 - End position of orf3a (Default = 26220)\n")
		print("9 - Start position of envelope protein (Default = 26245)\n")
		print("10 - End position of envelope protein (Default  = 26472)\n")
		print("11 - Start position of membrane glycoprotein (Default = 26523)\n")
		print("12 - End position of membrane glycoprotein (Default = 27191)\n")
		print("13 - Start position of orf6 (Default = 27202)\n")
		print("14 - End position of orf6 (Default = 27387)\n")
		print("15 - Start position of orf7a (Default = 27394)\n")
		print("16 - End position of orf7a (Default = 27759)\n")
		print("17 - Start position of orf7b (Default = 27756)\n")
		print("18 - End position of orf7b (Default = 27887)\n")
		print("19 - Start position of orf8 (Default = 27894)\n")
		print("20 - End position of orf8 (Default = 28259)\n")
		print("21 - Start position of nucleocapsid phosphoprotein (Default = 28274)\n")
		print("22 - End position of nucleocapsid phosphoprotein (Default = 29533)\n")
		print("23 - Start position of orf10 (Default = 29558)\n")
		print("24 - End position of orf10 (Default = 29674)\n")
		print("25 - Start position of vaccine (Default = 21508)\n")

		region_to_change = int(input())

		print("\nWhat is the new position?: \n")


		try:
			new_position = int(input())
			regions[region_to_change] = new_position
			print("\nChanged value:\n ------------\n " + str(regions[region_to_change]) + "\n")
		except:
			print("You must enter an integer value.")


	elif (user_action == 3):
		print("Start (index) of the vaccine: ?")
		vaccine_i = int(input())

	elif (user_action == 4):
		viralDNA = viralDNA(regions)
		virus = virus.replace("T", "U")
		vaccine = vaccine.replace("Ψ", "U")

		"""
		for i in range(len(virus)-len(vaccine)):
		  mm = virus[i:i+len(vaccine)]
		  mr = sum([c1 == c2 for c1,c2 in zip(mm, vaccine)]) / len(vaccine)
		  if mr > 0.5:
		    print(i, mr)
		exit(0)
		"""

		# vaccine starts at 21508 with a 67% match
		# spike protein starts at 21562
		vvirus = virus[vaccine_i:vaccine_i+len(vaccine)]

		print("\n Vaccine in the Virus Gene Seqence \n ------------------------ \n")
		print(vvirus)

		print("\n Vaccine Gene Seqence \n ------------------------ \n")
		print(vaccine)
		#viralDNA = viralDNA(untranslated_region_f = untranslated_region_f, orf1a_i = orf1a_i, orf1a_f = orf1a_f, orf1b_i = orf1b_i,orf1b_f = orf1b_f, spike_gp_i = spike_gp_i, spike_gp_f = spike_gp_f, orf3a_i = orf3a_i, orf3a_f = orf3a_f, envelope_p_i = envelope_p_i, envelope_p_f = envelope_p_f, membrane_gp_i = membrane_gp_i, membrane_gp_f = membrane_gp_f, orf6_i = orf6_i, orf6_f = orf6_f, orf7a_i = orf7a_i, orf7a_f = orf7a_f, orf7b_i = orf7b_i, orf7b_f = orf7b_f, orf8_i = orf8_i, orf8_f = orf8_f, n_p_i = n_p_i, n_p_f = n_p_f, orf10_i = orf10_i, orf10_f = orf10_f, new_region_i = 0, new_region_f = 10)


