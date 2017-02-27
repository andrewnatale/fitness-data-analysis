import os

# parse existing state lists
blacklists = {}
for listfile in os.listdir(os.getcwd()):
    if listfile.endswith(".lst"):
        # get structure
        structure_id = listfile[0:16]
        print structure_id
        if structure_id in blacklists:
            pass
        else:
            blacklists[structure_id] = []
        # parse list to find backrub ids
        with open(listfile, 'r') as state_list:
            reader = state_list.readlines()
        reader = [i.strip('\n') for i in reader]
        for entry in reader:
            brub_id = entry.split('/')[-1].split('_')[2][7:]
            print brub_id
            blacklists[structure_id].append(brub_id)

for structure_id in blacklists:
    with open('%s_blacklist' % structure_id, 'w') as outfile:
        for brub_id in sorted(blacklists[structure_id]):
            outfile.write('%s ' % brub_id)
