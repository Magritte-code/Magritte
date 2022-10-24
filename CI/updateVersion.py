#using regex to extract version number
import os
import re
import sys

# sys.argv[1] call with commit name
# if starts with the following: "major" for major, "minor" for minor and "patch" for patch (case insensitive)
# do the version update, otherwise just increase the patch version

this_dir = os.path.dirname(os.path.realpath(__file__))
with open(os.path.join(this_dir, "../version.txt"), "r+") as file:
    # Extract the lines containing the project description
    filecontents=file.read()
    # line = re.findall('project.*\([^\)]*\)', filecontents)[0]
    # print(line)
    # Extract the version number form those lines
    semver=list(re.finditer('(?P<major>\d*)\.(?P<minor>\d*)\.(?P<patch>\d*)', filecontents))

    versiondict=semver[0].groupdict()
    major=int(versiondict["major"])
    minor=int(versiondict["minor"])
    patch=int(versiondict["patch"])

    msg=sys.argv[1]#commit msg

    #matching word exactly (no exact letters directly afterwards, then any number of characters)
    #ofcourse, only the word itself may also be matched; removed end- and start-of-line tokens, as the commit messages can get too long
    major_reg=re.compile(r'.*\(increase version major\).*', re.IGNORECASE)
    minor_reg=re.compile(r'.*\(increase version minor\).*', re.IGNORECASE)
    patch_reg=re.compile(r'.*\(increase version patch\).*', re.IGNORECASE)
    #python 3.10 would enable us to use fancy switch statements; for now, we use if-else statements
    if re.match(major_reg, msg):
        #increment major version
        major=major+1
        minor=0
        patch=0
    elif re.match(minor_reg, msg):
        minor=minor+1
        patch=0
        #increment minor version
    elif re.match(patch_reg, msg):
        patch=patch+1
        #increment patch version
    else:
        patch=patch+1
        #by default, increment patch version

    # print([major, minor, patch])
    replaced_line=re.sub('(?P<major>|\d*)\.(?P<minor>\d*)\.(?P<patch>\d*)', str(major)+"."+str(minor)+"."+str(patch), filecontents)
    # print(replaced_line)

    #and finally replace the line in the file with the new line
    # line_regex=re.escape(line)
    # replaced_file=re.sub(line_regex, replaced_line, filecontents)
    # print("here")
    # print(replaced_file)

    file.close()

    writefile=open(os.path.join(this_dir, "../version.txt"), "w")

    writefile.write(replaced_line)
    writefile.close()

    #also output the new version
    sys.stdout.write(str(major)+"."+str(minor)+"."+str(patch))
