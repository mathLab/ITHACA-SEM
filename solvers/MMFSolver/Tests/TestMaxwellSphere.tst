<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>MMF Maxwell solver, DG, P=4</description>
    <executable>MMFSolver</executable>
    <parameters>MMFMaxwellSphere.xml</parameters>
    <files>
        <file description="Session File">MMFMaxwellSphere.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="F1" tolerance="1e-5">0.00208066</value>
            <value variable="F2" tolerance="1e-5">0.00208173</value>
            <value variable="F3" tolerance="1e-5">0.00180384</value>            
                        
        </metric>
        <metric type="Linf" id="2">
             <value variable="F1" tolerance="1e-5">0.00998094</value>
             <value variable="F2" tolerance="1e-5">0.00678984</value> 
             <value variable="F3" tolerance="1e-5">0.00225043</value>                                              
        </metric>
    </metrics>    
</test>


