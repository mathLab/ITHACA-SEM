<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Testing 3D homogeneous 2D derivatives</description>
    <executable>Deriv3DHomo1D</executable>
    <parameters>derivatives3Dhomo2D.xml</parameters>
    <files>
        <file description="Session File">derivatives3Dhomo2D.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-7">1.75249e-15</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-6">3.33067e-16</value>
        </metric>
    </metrics>
</test>


