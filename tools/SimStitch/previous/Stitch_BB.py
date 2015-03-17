import shutil
import tempfile
import os
import optparse
import traceback
import sys
import csv

from subprocess import check_call
from galaxy.jobs import JobDestination

class MockTool(object):

    def __init__(self, tool_dir):
        self.id = "client_test"
        self.version = "1.0"
        self.tool_dir = tool_dir

def main():
    """ Exercises a running lwr server application with the lwr client. """
    parser = optparse.OptionParser()
    parser.add_option('--outfile_log',dest='outfile_log',default='log.txt')
    parser.add_option('--outfile_stitch',dest='outfile_s',default='stitch_out.html')
    parser.add_option('--outdir_stitch', dest='outdir_s', default=None)
    parser.add_option('--outfile_peaks', dest='outfile_p', default = 'peaks_out.html')
    parser.add_option('--outdir_peaks',dest='outdir_p',default=None)
    parser.add_option('--infile', dest='infile',default=None)
    parser.add_option('--html_in', dest='html_in', default=None)
    parser.add_option('--html_in_dir', dest='html_in_dir', default=None)
    parser.add_option('--message_file', dest='message_file', default='messages.txt')
    parser.add_option('--noise_file', dest='noise_file', default='noise_level.txt')
    parser.add_option('--align_file', dest='align_file', default='alignment.txt')

	#So many Stitch parameters!
    parser.add_option('--cali', dest='txtFileCalibrants', default='')
    parser.add_option('--noisefilt', dest='noisefilt',default=1 )
    parser.add_option( '--incnoise', dest='incnoise',default=1)
    parser.add_option( '--minsnr', dest='minsnr',default=10)
    parser.add_option( '--remknown', dest='remknown',default=1)
    parser.add_option( '--knownnoise', dest='knownnoise',default="")
    parser.add_option( '--calon', dest='calon',default=1)	
    parser.add_option( '--calblank', dest='calblank',default=1)
    parser.add_option( '--calmode', dest='calmode',default=1)
    parser.add_option( '--calweighted', dest='calweighted',default=1)
    parser.add_option( '--calminrange', dest='calminrange',default=50)
    parser.add_option( '--calmaxpkd', dest='calmaxpkd',default=6.5)
    parser.add_option( '--calminsnr', dest='calminsnr',default=10)
    parser.add_option( '--calacm', dest='calacm',default=5e-5)
    parser.add_option( '--calacc', dest='calacc',default=0.002)
    parser.add_option( '--calacp', dest='calacp',default=0)
    parser.add_option( '--mzalign', dest='mzalign',default=1)
    parser.add_option( '--intcorr', dest='intcorr',default=0)
    parser.add_option( '--alignmin', dest='alignmin',default=20)
    parser.add_option( '--alignminsnr', dest='alignminsnr',default=6.5)
    parser.add_option( '--alignmaxpkd', dest='alignmaxpkd',default=1.5)
    parser.add_option( '--nullmzmin', dest='nullmzmin',default=0)
    parser.add_option( '--nullmzmax', dest='nullmzmax',default=0)
    parser.add_option( '--nullstart', dest='nullstart',default=15)
    parser.add_option( '--nullend', dest='nullend',default=15)
    parser.add_option( '--nullbnd', dest='nullbnd',default="70,2000")

    (options, args) = parser.parse_args()

    try:
        temp_directory = tempfile.mkdtemp()
        temp_work_dir = os.path.join(temp_directory, "w")

	for dir in [ temp_work_dir]:
            os.makedirs(dir)

        temp_input_path = os.path.join(temp_directory, "input.csv")

    except:
	error_file = open(options.outfile_s, 'w')
	error_file.write("Problem with creating temp files!")
	error_file.close()
	return


    try:

	with open(temp_input_path, "wb") as temp_input_file:
		csv_writer = csv.writer(temp_input_file,delimiter=',',quotechar=' ',quoting = csv.QUOTE_MINIMAL)
		csv_writer.writerow([2,0]) #THIS IS NOT ACTUALLY SET IN THE XML FILE. USERAWONLY_ON IS  NOT USUALLY ALTERED.
		csv_writer.writerow([3,options.noisefilt])
		csv_writer.writerow([4,options.incnoise])
		csv_writer.writerow([5,options.minsnr])
		csv_writer.writerow([6,options.remknown])
		csv_writer.writerow([7]+[options.knownnoise])	
		csv_writer.writerow([8,options.calon])
		csv_writer.writerow([9,options.calblank])
		csv_writer.writerow([10,options.calmode])
		csv_writer.writerow([11,options.calweighted])
		csv_writer.writerow([12,options.calminrange])
		csv_writer.writerow([13,options.calmaxpkd])
		csv_writer.writerow([14,options.calminsnr])
		csv_writer.writerow([15,options.calacm])
		csv_writer.writerow([16,options.calacc])
		csv_writer.writerow([17,options.calacp])
		csv_writer.writerow([18,options.mzalign])
		csv_writer.writerow([19,options.intcorr])
		csv_writer.writerow([20,options.alignmin])
		csv_writer.writerow([21,options.alignminsnr])
		csv_writer.writerow([22,options.alignmaxpkd])
		csv_writer.writerow([23,options.nullmzmin])
		csv_writer.writerow([24,options.nullmzmax])
		csv_writer.writerow([25,options.nullstart])
		csv_writer.writerow([26,options.nullend])
		csv_writer.writerow([27]+[options.nullbnd])
		csv_writer.writerow([28,0,0])
		#csv_writer.writerow([28,1,3]) # THIS IS NOT SET IN XML FILE. USED FOR MSMS TO SPECIFY WINDOW RANGE TO READ (user should use 0,0 or leave empty for full range to be used)
		csv_writer.writerow([29,options.txtFileCalibrants])

  	
	#return ##############################################################Just testing CSV file creation!
    	


    except:
	error_file = open(options.outfile_s, 'w')
	error_file.write("Problem with writing to temp files!")
	error_file.close()
	return

    command_inner =  "cd('../../../../tools/simstitch');Stitch_BB('{0}','{1}','{2}','{3}','{4}','{5}','{6}','{7}','{8}', '{9}', '{10}');exit;".format(options.infile,options.html_in, options.html_in_dir,temp_input_path, options.outfile_s,options.outdir_s, options.outfile_p, options.outdir_p, options.message_file, options.noise_file, options.align_file) 
    command_outer = "matlab -logfile \"{0}\" -r \"{1}\"".format(options.outfile_log, command_inner)	

    try:

	check_call([command_outer],shell=True)

    except:
	error_file = open(options.outfile_s, 'w')
	error_file.write(command_outer)
	error_file.close()
	return


    finally:
        shutil.rmtree(temp_directory)


if __name__ == "__main__":
    main()

