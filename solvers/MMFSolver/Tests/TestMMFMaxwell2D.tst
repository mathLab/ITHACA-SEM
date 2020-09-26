<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>MMF Maxwell solver, DG</description>
    <executable>MMFSolver</executable>
    <parameters>TestMMFMaxwell2D.xml</parameters>
    <files>
        <file description="Session File">TestMMFMaxwell2D.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="F1" tolerance="1e-5">5.31256e-05</value>
            <value variable="F2" tolerance="1e-5">5.2352e-05</value>
            <value variable="F3" tolerance="1e-5">4.0524e-05</value>            
                        
        </metric>
        <metric type="Linf" id="2">
             <value variable="F1" tolerance="1e-5">0.000315664</value>
             <value variable="F2" tolerance="1e-5">0.00050231</value> 
             <value variable="F3" tolerance="1e-5">0.000218756</value>                                              
        </metric>
    </metrics>    
</test>


