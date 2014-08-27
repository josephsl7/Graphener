'''
Created on Aug 26, 2014

@author: eswens13
'''

import Enumerator, Structs2Poscar, RunVasp, JobManager
import time
from sched import scheduler

if __name__ == '__main__':
    atomList = ['W','Re','Ni']
    
    enumerator = Enumerator.Enumerator(atomList)
    print "\nEnumerating symmetrically unique structures. . .\n"
    enumerator.enumerate()
    
    # It might be a good idea to have an 'Extractor' class once we figure out how to generate
    # the list of structures.
    enumerator.extract(['1','2','3'])
    
    print "\nConverting outputs to VASP inputs. . .\n"
    toPoscar = Structs2Poscar.Structs2Poscar(atomList)
    toPoscar.convertOutputsToPoscar()
    
    runner = RunVasp.RunVasp(atomList)
    print "\nPreparing directories for VASP. . .\n"
    runner.prepareForVasp()
    
    manager = JobManager.JobManager(atomList)
    
    s = scheduler(time.time, time.sleep)
    
    print "\nStarting low-precision ionic relaxation. . .\n"
    runner.run(1)
    
    finished = False
    start_time = time.time()
    event_time = start_time
    max_time = start_time + (8 * 60 * 60)
    while not finished and event_time <= max_time:
        event_time += 120
        s.enterabs(event_time, 1, manager.reportFinshed, ([runner.getCurrentJobIds()]))
        s.run()
        finished = manager.reportFinshed(runner.getCurrentJobIds())
    
    manager.reportLowStats()
    
    print "\nStarting normal-precision ionic relaxation. . .\n"
    runner.run(2)
    
    finished = False
    start_time = time.time()
    event_time = start_time
    max_time = start_time + (8 * 60 * 60)
    while not finished and event_time <= max_time:
        event_time += 120
        s.enterabs(event_time, 1, manager.reportFinshed, ([runner.getCurrentJobIds()]))
        s.run()
        finished = manager.reportFinshed(runner.getCurrentJobIds())
    
    manager.reportNormalStats()

    print "\nStarting DOS run. . .\n"
    runner.run(3)
    
    finished = False
    start_time = time.time()
    event_time = start_time
    max_time = start_time + (8 * 60 * 60)
    while not finished and event_time <= max_time:
        event_time += 120
        s.enterabs(event_time, 1, manager.reportFinshed, ([runner.getCurrentJobIds()]))
        s.run()
        finished = manager.reportFinshed(runner.getCurrentJobIds())
    
    manager.reportDOSStats()
    
        
    
    
    