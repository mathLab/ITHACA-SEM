<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Testing 3D homogeneous 2D derivatives</description>
    <executable>Deriv3DHomo2D</executable>
    <parameters>Deriv3D_Homo2D.xml</parameters>
    <files>
        <file description="Session File">Deriv3D_Homo2D.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value variable="dudx" tolerance="1e-7">3.05768e-14</value>
            <value variable="dvdy" tolerance="1e-7">0</value>
            <value variable="dwdz" tolerance="1e-7">0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="dudx" tolerance="1e-7">1.42109e-14</value>
            <value variable="dvdy" tolerance="1e-7">1.11022e-16</value>
            <value variable="dwdz" tolerance="1e-7">1.11022e-16</value>
        </metric>
    </metrics>
</test>


