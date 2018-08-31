<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 3D Helmholtz with cylindrical periodicity P=3 </description>
    <executable>ADRSolver</executable>
    <parameters>RotPerBcs3D.xml</parameters>
    <processes>4</processes>
    <files>
        <file description="Session File">RotPerBcs3D.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-10"> 1.42119e-07 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-10"> 4.06744e-05 </value>
        </metric>
    </metrics>
</test>
