infile = open('environment.yml')
outfile = open('environment_nobuilds.yml', 'w')
for line in infile:
    if '==' in line:
        outline = line.split('==')[0] + '\n'
    elif '=' in line:
        outline = line.split('=')[0] + '\n'
    else:
        outline = line
    outfile.write(outline)
infile.close()
outfile.close()