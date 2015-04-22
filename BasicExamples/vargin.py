#!/usr/bin/env python

#####
#vargin.py   
#Tue Apr 21 22:52:27 EDT 2015
#Ryan Melvin
#####
#This function take a input and output file path and returns them properly labeled.
#An example of getting variable input arguments. The following is apparently a standard C input parser.
#Also contains examples of
##Defining a python function
##Error catching
##For loop
##If statements
##String concatenation
##Changing behavior based on if the script is called by the user or something else.

#####
#Example call to python script
#python vargin.py -i ~/Desktop/input.dcd -o ~/Desktop.output.dat

#Successful output of this script
#Input file is  /Users/melvrl13/Desktop/input.dcd
#Ouput file is  /Users/melvrl13/Desktop.output.dat
#####
#
#####
#Procedures to Execute
#####

#The sys module from the python standard library allows access to command line arguments via sys.argv.
#The python standard library also has a module, getopt, for parsing input options based on standard identifies like '-', '--', '=' and ':'.
import sys, getopt

#Define a function 

def main(argv):
    #Initialize variables
    inputfile=''
    outputfile=''
    #Error catch for incorrect inputs
    try: #notice the indentation. That's how python knows where the "try" beings and ends
        #First, python will execute the code between "try:" and "except"
        #The following gets the options input by the user saved in a list "opts". 
        #In this example, there would be a option h with no value, a option i with a value required (indicated by ":") and a option o with a value required.
        #Opts will be of the form [(-option1, 'value1'),(-option2, 'value2'),....]. Args is anything input without an option preceding it.
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
        #If this code block executes successfully, the "except" block is skipped. Otherwise, an attempt is made to execute it.
    except getopt.GetoptError: 
        #Only an error from getopt will successfully execute this block. 
        #If another kind of error occurs, there will be an unhandled exception reported to the user.
        print 'vargin.py -i <inputfile> -o <ouputfile>'
        sys.exit(2) #The 2 here tells UNIX there was a command line syntax error.
        #Now we need to tell python what to do with each option
    for opt, arg in opts: #indenting is how you tell the script the loop has ended
        if opt == '-h':
                #This was our "help" option
                print 'test.py -i <inputfile> -o <outputfile>'
                sys.exit() #empty brackets indicate exiting normally
        elif opt in ("-i", "--ifile"): #The user may input a shorthand '-i' or longform '--ifile'
                inputfile = arg #assign to "inputfile" the value associated with the option "-i"
        elif opt in ("-o", "--ofile"): # I think you read this "if opt is a member of the set..."
                outputfile = arg
    print 'Input file is ', inputfile #Notice the indentation. This line is NOT in the for loop or if statement.
    print 'Ouput file is ', outputfile

#Checking who's in charge
#The "main" function is the one initiated by the user. So if this script is imported from another module, then "main" would be THAT module's name instead of this script. So we only ask "sys" what the input arguments were if the user was the one who called this function. In python, this type of checking lets us run certain parts when the code is called directly and other parts when it's called as part of something larger. 
if __name__ == "__main__":
    main(sys.argv[1:])


