<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Helmholtz solver in 2D domain</description>
    <executable python="true"> Helmholtz2D.py </executable>
    <parameters>Helmholtz2D_P7.xml</parameters>
    <files>
        <file description="Session File">../../MultiRegions/Tests/Helmholtz2D_P7.xml</file>
    </files>
    <metrics>
        <metric type="Linf" id="1">
            <value tolerance="1e-7">6.120209e-05</value>
        </metric>
    </metrics>
</test>
