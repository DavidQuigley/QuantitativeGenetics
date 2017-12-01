import os
import sys
import datetime
from datetime import date
import subprocess
import time
from optparse import OptionParser

class job:
    """Individual unit of work"""
    def __init__(self, sample_id, job_id,  job_dependencies, job_name, shell_command, is_multi_core):
        self.job_id = job_id
        self.dependencies = job_dependencies
        self.job_name = job_name
        self.shell_command = shell_command
        self.sample_id = sample_id
        self.is_multi_core = is_multi_core
    
    def intervals(self):
        if "INTERVAL" in self.parameters:
            return self.parameters["INTERVAL"]
        else:
            return []
    
    def write_command(self, fn_log):
        cmd = self.shell_command
        if fn_log != "":
            cmd =  'TIMESTAMP=$(date +"%Y-%m-%d %H:%M:%S")\n'
            cmd += 'echo "${TIMESTAMP} started job ' + self.job_name + ' id ' + str(self.job_id) + ' for ' + self.sample_id + '" >> ' + fn_log + "\n"
            cmd += self.shell_command
            cmd += 'TIMESTAMP=$(date +"%Y-%m-%d %H:%M:%S")\n'
            cmd += 'echo "${TIMESTAMP} completed job ' + self.job_name + ' id ' + str(self.job_id) + ' for ' + self.sample_id + '" >> ' + fn_log + "\n"
        
        return cmd

class jobDefinition:
    """
    stores shell command expansion and required parameters for a type of job
    One instance of this class is created for each kind of job read from the 
    job definition file. 
    Jobs that will request more than one core (multithreaded jobs) should be 
    created with job_multicore == True
    """
    def __init__(self, job_name, shell_command, required_parameters, job_multicore ):
        self.job_name = job_name
        self.shell_command = shell_command
        if not isinstance(required_parameters, list):
            print( "ERROR: required_parameters parameter must be a list of parameter names" )
            sys.exit(1)
        self.required_parameters = required_parameters  # list
        self.is_multi_core = job_multicore
    
    def expand_shell_command( self, parameters_passed ):
        """
        Instantiates a job definition with specific parameters
        Checks dict of passed parameters against known required parameters
        """
        if not isinstance(parameters_passed, dict):
            print( "ERROR: parameters_passed parameter must be a dict of parameter name->value" )
            sys.exit(1)
            
        for parameter in self.required_parameters:
            if not parameter in parameters_passed:
                print( "ERROR: required parameter " + parameter + " not passed to job " + self.job_name)
                sys.exit(1)
        
        # simple shell variable replacement; definitions can come either from 
        # global parameter list or individual parameters passed for a given job
        individual_command = self.shell_command
        for parameter in parameters_passed:
            individual_command = individual_command.replace( "${" + parameter + "}", parameters_passed[parameter] )
        
        return individual_command

class jobStub:
    """
    a jobStub is used for jobs that are used to define a workflow. 
    Each workflow needs a constant array of job stubs that can be re-used every time 
    that workflow is called by a execution. Stubs are turned into fully-realized jobs when instantiated
    """
    def __init__( self, job_id, dependency_ids, job_name, param_n2v ):
        self.job_name = job_name
        self.job_id = job_id
        self.dependency_ids = dependency_ids
        self.parameters = param_n2v

class workflowDefinition:
    """
    stores an array of job stubs. 
    """
    
    def __init__(self, workflow_name, parameters ):
        """
        workflow_name must matches a job defined with a jobDefinition, an internal 
        ID, an internal dependency list of zero or more items that is internally 
        consistent, and a list of parameters of the form key=value. All of the 
        parameters required by the job (as defined in the jobDefinition) should 
        be present as a key, and the values may include prefixes or suffixes 
        (e.g. instead of FN_OUT=${FN_OUT} individual jobs might use 
        FN_OUT=${FN_OUT}_sorted.bam to make the whole workflow work together.
        """
        self.name = workflow_name
        self.parameter_list = parameters
        self.job_stubs = []
    
    def add_job( self, job_id, upstream_job_ids, job_name, parameters, job_definitions, line_no ):
        
        """
        adds a new job to the internal array of job stubs
        check that job_id and upstream job ids are integers, comma-sep. integers respectively
        """
        try:
           job_id = int( job_id )
        except ValueError as err:
            print( "ERROR: workflow template error on line " + str(line_no) + ": job_id must be integer")
            sys.exit(1)
        try: 
            upstream_job_ids = [ int(x) for x in upstream_job_ids.split(",") ]
        except ValueError as err:
            print( "ERROR job dependencies on workflow line " + str(line_no) + \
                    " must be comma-separated list of integers")
            sys.exit(1)
        
        # upstream_job_ids equal to 0 means no dependencies. already converted to integer.
        if len(upstream_job_ids)==1 and upstream_job_ids[0]==0:
            upstream_job_ids=[]
                                
        # check that job_id is not already present
        for stub in self.job_stubs:
            if job_id == stub.job_id:
                print( "ERROR: cannot allow duplicate job ID in workflow with name " + job_name )
                sys.exit(1)
        
        # check that dependencies are already defined
        for upstream_job_id in upstream_job_ids:
            seen = False
            for stub in self.job_stubs:
                if upstream_job_id == stub.job_id:
                    seen = True
            if not seen:
                print ("ERROR: dependency on job " + str(upstream_job_id) + " in job " + \
                        str(job_id) + " in job " + job_name + \
                        " on workflow line " + str(line_no) + " declared before it was defined")
                sys.exit(1)
        
        # check that this job has a definition 
        if not job_name in job_definitions:
            print( "ERROR: workflow refers to " + undefined job " + job_name)
            sys.exit(1)
        # check that required parameters will be filled by the stub
        param_n2v = {}
        
        for parameter in [x.strip() for x in parameters.split(",")]:
            try:
                key, value = parameter.split("=")
                param_n2v[key] = value
            except ValueError as err:
                print( "ERROR: job definition syntax error " + str(line_number) + ", missing comma in key-value pair")
                sys.exit(1)
        
        for required_parameter in job_definitions[job_name].required_parameters:
            if not required_parameter in param_n2v:
                print( "ERROR: job definition for " + job_name + \
                  " requires parameter " + required_parameter + " not " + \
                  " present in workflow definition for " + self.name )
                sys.exit(1)
        stub = jobStub( job_id, upstream_job_ids, job_name, param_n2v )
        self.job_stubs.append( stub )

        
class workflowManager:
    
    def __init__(self, verbose, memory_per_core, memory_monolithic, n_cores ):        
        self.job_definitions = {}
        self.workflow_definitions = {}
        self.jobs = []
        self.environment = ""
        self._workflow_id2highest_external_id = {}   # for each workflow, highest external id
        self._verbose = verbose
        self.memory_monolithic = memory_monolithic # max mem, single core 
        self.memory_per_core = memory_per_core     # max mem per core, multicore 
        self.n_cores = n_cores                     # number of cores to require for multicore
        self.dir_scripts = ""
        self.dir_logs = ""
    
    def parse_job_definitions( self, fn ):
        """
        parse a file that defines individual jobs
        definitions starts with #job JOBNAME {SINGLE_CORE,MULTI_CORE}, VARIABLES
        where only one of {SINGLE_CORE,MULTI_CORE} is passed and VARIABLES 
        is a comma-delimited list of required variables. SAMPLE_ID is always 
        a required variable.
        
        Example: 
        #job hello SINGLE_CORE SAMPLE_ID
        echo "Hello world" > /tmp/hello.txt
        #endjob
        
        Parameter substitution for jobs can come either from environment variables  
        common to all parts of the workflow or parameters are passed in from a 
        particular workflow instantiation of a workflow
        """
        
        f = open( fn )
        reading_shell_command=False
        shell_command = ""  
        line_no=0      
        for line in f:
            line_no += 1
            if line[0] == "#":
                a = line[1:].rstrip('\r\n').split(' ')
                if a[0]=="endjob":
                    self.job_definitions[ job_name ] = jobDefinition( job_name, shell_command, job_parameters, job_multicore )
                    shell_command = ""
                    reading_shell_command = False
                
                elif a[0] == "job":
                    job_name = a[1]
                    job_core = a[2]
                    if job_core!="MULTI_CORE" and job_core!="SINGLE_CORE":
                        print("ERROR: second element of job definition line " \
                               + str(line_no) + " is " + job_core + \
                               " but must be one of {MULTI_CORE,SINGLE_CORE}")
                        sys.exit(1)
                    job_multicore = job_core=="MULTI_CORE"
                    job_parameters = a[3].split(',')
                    reading_shell_command = True
            
            elif reading_shell_command:
                line_to_add = line.rstrip('\r\n')
                if line_to_add == "":
                    line_to_add = "\n"
                shell_command += line_to_add + '\n'
        
        f.close()
        if self._verbose:
            print "MARLOWE: Parsed job types " + ", ".join( self.job_definitions.keys() )
    

    def parse_workflow_definitions( self, fn ):
        """
        A workflow is a series of jobs that together accomplish a bigger goal, 
        such as "align a FASTQ file, deduplicate, and index".
        A workflow definition file contains one or more workflows.
        The workflow definition file has the format
        
        #workflow  workflow_name PARAM_KEY_,PARAM_KEY_2...PARAM_KEY_N
        jobid	dep,dep,dep	job_name	PARAM_KEY_1=PARAM_VAL_1,...
        #endworkflow
        
        The first line defining the workflow has 3 tab-separated elements: #workflow,
        the name of the workflow, and a comma-separated list of parameter keys 
        that must be passed in by a particular workflow. 
        subsequent lines declare individual jobs that have an internal job 
        ID (i.e. one that is only to manage dependencies within the workflow).
        The last line #endworkflow ends the workflow definition. Definitions 
        can be separated by any number of lines.
        """
        f = open( fn )
        reading=False
        line_no=0
        
        for line in f:
            line_no += 1
            if line[0] == "#":
                a = line[1:].rstrip('\r\n').split('\t')
                if a[0]=="endworkflow":
                    self.workflow_definitions[ workflow_name ] = workflow
                    reading=False
                
                elif a[0]=="workflow":
                    if len(a) != 3:
                        print("ERROR: workflow definition on line " + \
                               str(line_no) + " does not have 3 tab-separated elements")
                        sys.exit(1)
                    js, workflow_name, parameter_list = a
                    parameter_list = [ x.strip() for x in parameter_list.split(",") ]
                    workflow = workflowDefinition( workflow_name, parameter_list )
                    reading=True # now reading individual jobs
                
            elif reading:
                if line.rstrip('\r\n') == "":
                    continue
                a = line.rstrip('\r\n').split('\t')
                if len(a) != 4:
                    print("ERROR: did not find four tab-separated elements for workflow " + workflow_name + " on line " + str(line_no) )
                    sys.exit(1)
                job_id, job_dependency_list, job_name, parameters_passed = ( a[0], a[1], a[2], a[3] )
                # sanity checking takes place in add_job()
                workflow.add_job( job_id, job_dependency_list, job_name, parameters_passed, self.job_definitions, line_no )
        f.close()
        if reading:
            print ("ERROR: workflow not closed at end of workflow definition file")
            sys.exit(1)
        
        if self._verbose:
            print "# MARLOWE parsed workflow types " + ", ".join( self.workflow_definitions.keys() )
    
    
    def parse_execution( self, fn ):
        """
        The execution is a particular call to instantiate one or more workflows.
        Individual execution items have an ID number, and can depend on other 
        previously defined execution calls. Every line in the execution definition 
        file has the form
        # MyID	dep_IDs	workflow key=value,key=value...
        dep_IDs is a comma-delimited list of ID numbers on which a given item 
        depends. If an item does not depend on any previously defined item, give  
        it a dep_IDs value of 0. The value of workflow much be defined in the 
        workflow definition file. Keys defined in the workflow will be passed on 
        to the workflow when it is instantiated.
        """
        f = open( fn )
        line_number=0
        execution_ids = []
        for line in f:
            line_number += 1 
            if line.rstrip('\r\n')=="" or line[0] == "#":
                pass
            else:
                a = line.rstrip('\r\n').split("\t")
                if len(a) != 4:
                    print("ERROR: did not find four tab-separated elements for job workflow on line " + str(line_number) )
                    sys.exit(1)
                execution_id, execution_dependency_list, workflow_name, parameters_passed = ( a[0], a[1], a[2], a[3] )
                if not workflow_name in self.workflow_definitions: 
                    print("ERROR: execution specifies undefined workflow named " + workflow_name + ", check workflow template file" )
                    sys.exit(1)
                try:
                    execution_id = int( execution_id )
                except ValueError as err:
                    print( "ERROR: execution identifier on line " + str(line_number) + " of execution file is not an integer")
                    sys.exit(1)
                try: 
                    upstream_execution_ids = [ int(x) for x in execution_dependency_list.split(",") ]
                except ValueError as err:
                    print( "ERROR: execution dependencies on line " + str(line_number) + " of execution file must be comma-separated list of integers")
                    sys.exit(1)
                if len(upstream_execution_ids)==1 and upstream_execution_ids[0]==0:
                    upstream_execution_ids=[]
                else:
                    for upstream_execution_id in upstream_execution_ids:
                        if not upstream_execution_id in execution_ids:
                            print("ERROR: dependency on execution id " + str(upstream_execution_id) + " in " + \
                                 str(execution_id) + " declared on line " + str(line_number) + " of execution file before it was defined")
                            sys.exit(1)
                # now know name of workflow, execution_id, [upstream_execution_ids], and parameters. No substitution here.
                workflow_param_n2v = {}
                for parameter in parameters_passed.split(','):
                    try:
                        parameter_name, parameter_value = parameter.split("=")
                    except ValueError:
                        print( "ERROR: parsing comma-delimited list of parameters with format name=value in job " + str(execution_id) )
                        sys.exit(1)
                    workflow_param_n2v[parameter_name.strip()] = parameter_value.strip()
                
                # we now have all of the information required to create the actual jobs
                # instantiating the jobs will, for each job in the workflow, assign 
                # an external job ID, external depndencies, write a shell command
                execution_ids.append(execution_id)
                self.instantiate_jobs( execution_id, upstream_execution_ids, workflow_name, workflow_param_n2v )
    
    def instantiate_jobs( self, workflow_id, upstream_workflow_ids, workflow_name, workflow_param_n2v ):
        
        workflow = self.workflow_definitions[ workflow_name ]
        # the workflow has a list of jobs. 
        # create a set of Job objects that contain external job IDs. To calc 
        # external ID, for each job increment job ID and all non-zero dependency
        # values by the number of jobs currently assigned. To account for 
        # dependencies between workflows, if this workflow has any dependencies, 
        # assign individual jobs where the dependency list equals zero to be 
        # the highest value from the upstream workflows
        
        self._workflow_id2highest_external_id[workflow_id]=0
        
        external_increment = len(self.jobs)
        for stub in workflow.job_stubs:
            job_id = stub.job_id + external_increment
            if job_id > self._workflow_id2highest_external_id[workflow_id]:
                self._workflow_id2highest_external_id[workflow_id]=job_id
            job_dependencies = []
            if len( stub.dependency_ids ) > 0:
                # has internal dependencies
                for id in stub.dependency_ids:
                    job_dependencies.append( id + external_increment )
            else:
                # no internal dependencies
                if len( upstream_workflow_ids ) > 0:
                    # but external dependencies
                    job_dependencies = []
                    for id in upstream_workflow_ids:
                        job_dependencies.append( self._workflow_id2highest_external_id[id] )
            job_name = stub.job_name
            
            # substitute in parameter values from the contents of workflow_param_n2b 
            
            stub_param_n2v = dict(stub.parameters)
            for key_stub in stub_param_n2v:
                for key_work in workflow_param_n2v:
                    target = "${" + key_work + "}"
                    replacement = workflow_param_n2v[key_work]
                    if target in stub_param_n2v[key_stub]:
                        stub_param_n2v[key_stub] = stub_param_n2v[key_stub].replace( target, replacement )
            
            if not "SAMPLE_ID" in stub_param_n2v:
                print("ERROR: required parameter SAMPLE_ID not defined in workflow " + workflow_name)
                sys.exit(1)
            sample_id = stub_param_n2v["SAMPLE_ID"]
            shell_command = self.job_definitions[ job_name ].expand_shell_command( stub_param_n2v )
            is_multi_core = self.job_definitions[ job_name ].is_multi_core
            newjob = job( sample_id, job_id, job_dependencies, job_name, shell_command, is_multi_core)
            self.jobs.append( newjob )
        
    def write_jobs_to_queue( self, workflow_id ):
        ### Job dependencies are handled through the hold_jid parameter of qsub
        
        fn_logfile=self.dir_logs + "/" + workflow_id + "_log.txt"
        f_log = open(fn_logfile, "w")
        f_log.write( "Began submitting jobs to queue\n")
        for job in self.jobs:
            job_identifier = "j" + str(job.job_id) + "_" + workflow_id 
            f_log.write( datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + " submitted " + job_identifier + " to execute " + job.job_name + "\n" )
            # create job script, push job script onto queue
            script = "#!/bin/bash\n"
            script += "# created " + str( datetime.datetime.now() ) + "\n"
            script += "MAX_MEM_PER_CORE=" + str(self.memory_per_core) + "G\n"
            script += "MAX_MEM_TOTAL=" + str(self.memory_monolithic) + "G\n"
            script += "N_CORES=" + str(self.n_cores) + "\n"
            script += self.environment
            if True:#job.is_multi_core:
                script += "#$ -pe orte " + str(self.n_cores) + '\n'
            mem_request = self.memory_monolithic
            if job.is_multi_core and (self.n_cores * self.memory_per_core) > mem_request:
                mem_request = self.n_cores * self.memory_per_core
            script += "#$ -l h_vmem=" + str(mem_request) + "G\n"
            script += "#$ -N " + job_identifier + "\n"
            if len( job.dependencies ) > 0:
                str_deps = ",".join( [ "j" + str(x) + "_" + workflow_id for x in job.dependencies ] )
                script += "#$ -hold_jid " + str_deps + '\n'
            script += "#$ -e " + self.dir_scripts + "/" + job_identifier + "_output.log\n"
            script += job.write_command( fn_logfile )
            path_to_write = self.dir_scripts + '/' + job_identifier + '.sh'
            fo = open( path_to_write, "w")
            fo.write( script )
            fo.close()
            try:
                subprocess.call( ["qsub", path_to_write] )
            except OSError as err:
                print "ERROR: qsub not found, is Sun Grid Engine installed?"
                sys.exit(1)
        
        f_log.close()
    
    def write_jobs_to_stdout( self, workflow_id ):
        """Print a dry run of jobs and the queue assignments to stdout"""
        
        fn_logfile = ""
        if self.dir_logs != "":
            fn_logfile=self.dir_logs + "/" + workflow_id + "_log.txt"            
        
        print( self.environment )
        print( "N_CORES=" + str(self.n_cores)  )
        print( "MAX_MEM_PER_CORE=" + str(self.memory_per_core) + "G" )
        print( "MAX_MEM_TOTAL=" + str(self.memory_monolithic) + "G" )
        for job in self.jobs:
            print( job.write_command( fn_logfile ) )
        print( "#------------------------------------------------" )
        print( "# jobs in queue" )
        print( "#------------------------------------------------" )
        for job in self.jobs:
            s = "# " + str(job.job_id)
            if len( job.dependencies ) > 0:
                s = s + " " + ",".join( [ str(x) for x in job.dependencies ] )
            else:
                s = s + " 0"
            print( "#" + s + " " + job.job_name + " " + job.sample_id )
            
    
    def read_environment( self, fn_environment ):
        """Read file containing environment variables to be printed verbatim"""
        f = open( fn_environment )
        self.environment = ''.join( f.readlines() )
        f.close()

#-------------------------------------------------------------------------------

default_execution=datetime.datetime.now().strftime( "%y%m%d_%H%M%S" )

parser = OptionParser()
parser.add_option("-t", "--fn_job_template", dest="fn_job_templates", \
                  help="REQUIRED file containing individual job templates", default="")
parser.add_option("-j", "--fn_workflow_templates", dest="fn_workflow_templates", \
                  help="REQUIRED file containing workflow templates", default="")
parser.add_option("-f", "--fn_execution", dest="fn_execution", \
                  help="REQUIRED file containing execution to submit to queue", default="")
parser.add_option("-e", "--fn_environment", dest="fn_environment", \
                  help="optional, file containing environment variables ", default="")                  
parser.add_option("-s", "--dir_scripts", dest="dir_scripts", \
                  help="REQUIRED directory to write scripts that will be executed on queue", \
                  default="")
parser.add_option("-p", "--execution_prefix", dest="execution_prefix", \
                  help="execution prefix, default to date/time string", \
                  default=default_execution)
parser.add_option("-l", "--dir_logs", dest="dir_logs", \
                  help="optional, directory to write log files; defaults to value of dir_scripts", \
                  default="")
parser.add_option("-g", "--memory_per_core", dest="mem_per_core", default="2", \
                  help="Gb of RAM to require per core for multithreaded jobs (e.g. BWA)")
parser.add_option("-n", "--n_cores", dest="n_cores", default="1", \
                  help="Number of cores to require for jobs that can run multithreaded")
parser.add_option("-m", "--memory_monolithic", dest="mem_monolithic", default="2", \
                  help="Gb of RAM to require for non-threadable jobs (e.g. GATK)")
parser.add_option("-r", "--dry_run", action="store_true", dest="dry_run", \
                  help="write execution script to stdout but do not run", \
                  default=False)
parser.add_option("-x", "--timestamp", action="store_true", dest="timestamp", \
                  help="emit timestamp", \
                  default=False)                  
parser.add_option("-v", "--verbose", action="store_true", dest="verbose", \
                  help="write progress and diagnostic messages to stdout", default=False)

(options, args) = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    sys.exit(1)

fn_job_templates = options.fn_job_templates
fn_workflow_templates = options.fn_workflow_templates
fn_execution = options.fn_execution
fn_environment = options.fn_environment
dir_logs = options.dir_logs
dir_scripts = options.dir_scripts
mem_per_core = options.mem_per_core
mem_monolithic = options.mem_monolithic
n_cores = options.n_cores
timestamp = options.timestamp

if timestamp:
    print( str( datetime.datetime.now() ) )
    sys.exit(0)
    
try:
    n_cores = int(n_cores)
except ValueError as err:
    print "value of --n_cores parameter must be an integer"
    sys.exit(1)
try:
    mem_per_core = int(mem_per_core)
except ValueError as err:
    print "value of --memory_per_core parameter must be an integer"
    sys.exit(1)
try:
    mem_monolithic = int(mem_monolithic)
except ValueError as err:
    print "value of --memory_monolithic parameter must be an integer"
    sys.exit(1)

dry_run = options.dry_run
execution_prefix = options.execution_prefix
verbose = options.verbose

if verbose:
    print "# MARLOWE reading job template file: " + fn_job_templates
    print "# MARLOWE reading workflow template file: " + fn_workflow_templates
    print "# MARLOWE reading execution file: " + fn_execution

if not os.access( fn_job_templates, os.R_OK ):
    print( "ERROR: Cannot read from job template file specified in --fn_job_template" )
    sys.exit(1)
if not os.access( fn_execution, os.R_OK ):
    print( "ERROR: Cannot read from execution file specified in --fn_execution" )
    sys.exit(1)
if not os.access( fn_workflow_templates, os.R_OK ):
    print( "ERROR: Cannot read from workflow template file specified in --fn_workflow_templates" )
    sys.exit(1)
if fn_environment != "" and not os.access( fn_environment, os.R_OK ):
    print( "ERROR: Cannot read from environment variable file specified in --fn_environment" )
    sys.exit(1)
    
manager = workflowManager( verbose, mem_per_core, mem_monolithic, n_cores )
manager.dir_scripts = dir_scripts
manager.dir_logs = dir_logs
manager.parse_job_definitions( fn_job_templates )
manager.parse_workflow_definitions( fn_workflow_templates )
manager.parse_execution( fn_execution )

if fn_environment != "":
    manager.read_environment( fn_environment )

if verbose:
    print "# MARLOWE writing scripts to: " + dir_scripts
    print "# MARLOWE writing logs to: " + dir_logs
    print "# MARLOWE read " + str(len(manager.jobs)) + " job(s) from execution"

if dry_run:
    manager.write_jobs_to_stdout(execution_prefix)
else:
    manager.dir_scripts = dir_scripts
    if not os.path.exists( manager.dir_scripts ):
        print( "ERROR: script output directory specified in --dir_scripts does not exist" )
        sys.exit(1)    
    if not os.access( manager.dir_scripts, os.W_OK ):
        print( "ERROR: cannot write to script output directory specified in --dir_scripts")
        sys.exit(1)
    
    manager.write_jobs_to_queue(execution_prefix)

