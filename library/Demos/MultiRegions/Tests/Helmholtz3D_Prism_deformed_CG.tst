<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Helmholtz 3D CG for deformed Prism</description>
    <executable>Helmholtz3D</executable>
    <parameters>helmholtz3D_prism_deformed.xml</parameters>
    <files>
        <file description="Session File">helmholtz3D_prism_deformed.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-9">1.61137e-06</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-9">8.98093e-06</value>
        </metric>
    </metrics>
</test>


