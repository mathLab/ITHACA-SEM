<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Helmholtz 3D HDG for Prism</description>
    <executable>HDGHelmholtz3D</executable>
    <parameters>--use-scotch Helmholtz3D_Prism.xml</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">Helmholtz3D_Prism.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-9">0.000205416</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-9">0.00322597</value>
        </metric>
    </metrics>
</test>
