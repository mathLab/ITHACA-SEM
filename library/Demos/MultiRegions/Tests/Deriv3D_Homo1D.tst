<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Testing 3D homogeneous 1D derivatives</description>
    <executable>Deriv3DHomo1D</executable>
    <parameters>Deriv3D_Homo1D.xml</parameters>
    <files>
        <file description="Session File">Deriv3D_Homo1D.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value variable="dudx" tolerance="1e-7">4.52931e-15</value>
            <value variable="dvdy" tolerance="1e-7">0.00406798</value>
            <value variable="dwdz" tolerance="1e-7">2.86856e-16</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="dudx" tolerance="1e-7">3.55271e-15</value>
            <value variable="dvdy" tolerance="1e-7">0.00464427</value>
            <value variable="dwdz" tolerance="1e-7">2.22045e-16</value>
        </metric>
    </metrics>
</test>


