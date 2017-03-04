<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>MMF Maxwell solver, DG, P=4</description>
    <executable>MMFSolver</executable>
    <parameters>TestMaxwellSphereh07.xml</parameters>
    <files>
        <file description="Session File">TestMaxwellSphereh07.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="F1" tolerance="1e-5">0.108817</value>
            <value variable="F2" tolerance="1e-5">0.129406</value>
            <value variable="F3" tolerance="1e-5">0.0580759</value>            
                        
        </metric>
        <metric type="Linf" id="2">
             <value variable="F1" tolerance="1e-5">0.146007</value>
             <value variable="F2" tolerance="1e-5">0.164459</value> 
             <value variable="F3" tolerance="1e-5">0.063555</value>                                              
        </metric>
    </metrics>    
</test>


