<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>3D Helmholtz with cylindrical periodicity, annulus, P=3</description>
    <executable>ADRSolver</executable>
    <parameters>RotPerBcs3D_Annulus.xml</parameters>
    <processes>4</processes>
    <files>
        <file description="Session File">RotPerBcs3D_Annulus.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-10">2.35947e-10</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-10">1.05215e-09</value>
        </metric>
    </metrics>
</test>
