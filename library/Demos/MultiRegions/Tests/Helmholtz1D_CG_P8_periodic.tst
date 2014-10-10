<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Helmholtz 1D CG with P=8, periodic BCs</description>
    <executable>Helmholtz1D</executable>
    <parameters>Helmholtz1D_P8_periodic.xml</parameters>
    <files>
        <file description="Session File">Helmholtz1D_P8_periodic.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-9">4.00482e-10</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-9">1.64741e-10</value>
        </metric>
    </metrics>
</test>


