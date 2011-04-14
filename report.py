#!/usr/bin/python

import sys,os,re,shutil

filelist = ["logarithmic_form.py","logarithmic_forms.py","singular_module.py",
                            "graded_module.py","logarithmic_derivations.py"]

def get_uid():
    #Hacky concurrency - need to set this up bf project
    uid_file = open(".uid","r")
    uid = int(uid_file.read())
    uid_file.close()
    uid_file = open(".uid","w")
    uid_file.write(str(uid+1))
    uid_file.close()
    return uid

def make_dir(uid):
    os.mkdir("reports/report_"+str(uid))

def get_dir(uid):
    return "reports/report_"+str(uid)

def copy_env(uid):
    for pyfile in filelist:
        shutil.copy(pyfile,os.path.join("reports","report_"+str(uid)))

def vars_from_divisor(d_string):
    var_list = []
    underscore_vars = re.findall("[a-z]_\d+",d_string)
    single_vars = re.findall("^[a-z][a-z]^[a-z_]",d_string)
    first_single = re.findall("[a-z][^a-z_]",d_string)
    last_single = re.findall("[^a-z_][a-z]",d_string)
    if len(first_single)>0 and len(first_single[0])>0:
        var_list.append(first_single[0][0])
    var_list.extend([ls[1] for ls in last_single])
    var_list.extend(underscore_vars)
    for u_v in underscore_vars:
        if u_v[0] in var_list:
            var_list.remove(u_v[0])
    var_list.extend([sv[1] for sv in single_vars])
    print var_list
    return list(set(var_list))

def create_divisor(d_string):
    var_list = vars_from_divisor(d_string)
    poly_str = "C.<"+str(",".join(var_list))+"> = PolynomialRing(QQ,"+str(len(var_list))+")\n"
    divisor_str = "divisor = "+d_string+"\n"
    return poly_str+divisor_str

def write_include(uid,num,divisor):
    include_file = open(os.path.join("reports","report_"+str(uid),"report"+str(i)+".tex"),"w")
    include_master = open("report_base.tex","r")
    include_string = include_master.read()
    include_master.close()
    include_string = include_string.replace("%Divisor%",create_divisor(divisor))
    include_string = include_string.replace("%DivisorName%",divisor)
    include_file.write(include_string)
    include_file.close()

def write_master(uid,d_strings):
    master_file = open(os.path.join("reports","report_"+str(uid),"master_report.tex"),"w")
    master_master = open("report_master.tex","r")
    master_string = master_master.read()
    master_master.close()
    inputs = []
    for i,d_string in enumerate(d_strings):
        input_file = os.path.join("report"+str(i)+".tex")
        inputs.append("\input{"+input_file+"}")
    inputs = "\n".join(inputs)
    master_string = master_string.replace("%Inputs%",inputs)
    master_file.write(master_string)
    master_file.close()

cwd = os.getcwd()
d_strings = sys.argv[1:]
uid = get_uid()
make_dir(uid)
write_master(uid,d_strings)
for i,d_string in enumerate(d_strings):
    write_include(uid,i,d_string)
copy_env(uid)
os.chdir(get_dir(uid))
os.system("pdflatex master_report.tex -no-shell-escape")
os.system("sage master_report.sage")
os.system("pdflatex master_report.tex -no-shell-escape")
os.system("pdflatex master_report.tex -no-shell-escape") # Biog
os.chdir(cwd)

