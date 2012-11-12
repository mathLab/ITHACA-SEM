<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Testing 3D homogeneous 1D derivatives</description>
    <executable>Deriv3DHomo1D</executable>
    <parameters>derivatives3Dhomo1D.xml</parameters>
    <files>
        <file description="Session File">derivatives3Dhomo1D.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-7">6.40541e-15</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-6">3.55271e-15</value>
        </metric>
    </metrics>
</test>


