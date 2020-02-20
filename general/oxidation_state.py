import sys

oxidation_state = {'Pd':[0,2,4],'O':[-2]} #oxidation state dictionary
composition     = {}

'Unique compositions'
for i in [1,2]:
	for j in [1,2]:
		name = 'Pd'+str(j)+'O'+str(i)
		#name = name.replace('1','')
		composition[name] = []

'oxidation states of each composition'
for comp in composition:
	for ox1 in oxidation_state[comp[0:2]]:
		for ox2 in oxidation_state[comp[-2:-1]]:
			composition[comp].append(int(ox1)*int(comp[2])+int(ox2)*int(comp[-1]))

'pure metal'
for element in oxidation_state:
	if element not in ['O','H']:
		composition[element] = oxidation_state[element]

print(composition)

print('Oxidation States = 0,+1+,2:')
'oxidation states in 0,+1,+2'
for comp in composition:
	for ox in composition[comp]:
		if ox in [0,1,2]:
			print(comp.replace('1',''), ox)
