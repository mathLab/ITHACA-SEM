<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Helmholtz 1D CG with P=8 and Robin BCs</description>
    <executable>Helmholtz1D</executable>
    <parameters>helmholtz1D_8modes_RBC.xml</parameters>
    <files>
        <file description="Session File">helmholtz1D_8modes_RBC.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-7">1.11262e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-6">1.28659e-05</value>
        </metric>
    </metrics>
</test>


