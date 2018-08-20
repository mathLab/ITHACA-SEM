<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Helmholtz 3D CG for deformed Prism</description>
    <executable>Helmholtz3D</executable>
    <parameters>Helmholtz3D_Prism_Deformed.xml</parameters>
    <files>
        <file description="Session File">Helmholtz3D_Prism_Deformed.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-9">8.71461e-06</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-9">7.59576e-05</value>
        </metric>
    </metrics>
</test>


