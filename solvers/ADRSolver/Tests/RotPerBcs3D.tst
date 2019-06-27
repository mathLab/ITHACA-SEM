<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>3D Helmholtz with cylindrical periodicity P=3</description>
    <executable>ADRSolver</executable>
    <parameters>RotPerBcs3D.xml</parameters>
    <files>
        <file description="Session File">RotPerBcs3D.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-10">1.56506e-11</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-10">5.01964e-09</value>
        </metric>
    </metrics>
</test>
