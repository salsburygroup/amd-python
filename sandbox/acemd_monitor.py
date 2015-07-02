#!/usr/env/ Python

#Ryan Melvin
#Thu Jul  2 16:38:56 EDT 2015

#Run acemd, monitor for termination and email user last 9 lines of log file

import subprocess
import smtplib
from email.mime.text import MIMEText
import argparse
import socket

# Initialize parser for user input
parser = argparse.ArgumentParser(description ='Run acemd, monitor for termination and email user last 9 lines of log file', add_help=False)
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-d', '--device', action='store', dest='device', help='Which GPU card to use',type=str,required=True)
inputs.add_argument('-c', '--config', action='store', dest='conf', help='configuration file', type=str, required=True)
inputs.add_argument('-l', '--log', action='store', dest='logfile', help='Log file. Assumes cwd if no other path specified',default='log.txt')

UserInput=parser.parse_args()

acemd_cmd = 'acemd --device ' + UserInput.device + ' ' + UserInput.conf + ' >&' + UserInput.logfile

process = subprocess.Popen(acemd_cmd, shell=True)
process.wait()

# Get Machine IP address for email subject
# Seems overdone but makes sure you don't get the linux default 127.0.0.1
check_ip = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
check_ip.connect(('google.com', 0))
ip_address = check_ip.getsockname()[0]

# Open the log file
fp = open(UserInput.logfile,'rb')

#Prepare last 9 lines for email 
msg = MIMEText("\n".join(fp.read().strip().split("\n")[-9:]))
fp.close() #Close the file

# Prepare message for smtp
sender = 'melvrl13@wfu.edu'
recipient = 'ryanlmelvin@gmail.com'

msg['Subject'] = UserInput.conf + ' on ' + ip_address + ' device ' + UserInput.device 
msg['From'] = sender
msg['To'] = recipient

# Credentials
username = 'melvrl13@wfu.edu'
password = 'thisis4th'

s = smtplib.SMTP('smtp.gmail.com:587')
s.starttls()
s.login(username,password)
s.sendmail(sender, recipient, msg.as_string())
s.quit()
