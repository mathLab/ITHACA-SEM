<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 3D Helmholtz with cylindrical periodicity P=4 </description>
    <executable>ADRSolver</executable>
    <parameters>RotPerBcs3D.xml</parameters>
    <processes>4</processes>
    <files>
        <file description="Session File">RotPerBcs3D.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12"> 3.8907e-09 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-10"> 1.39131e-06 </value>
        </metric>
    </metrics>
</test>
