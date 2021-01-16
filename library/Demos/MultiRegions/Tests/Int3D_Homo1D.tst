<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Testing 3D homogeneous 1D integral</description>
    <executable>Int3DHomo1D</executable>
    <parameters>Int3D_Homo1D.xml</parameters>
    <files>
        <file description="Session File">Int3D_Homo1D.xml</file>
    </files>

    <metrics>
        <metric type="Linf" id="2">
            <value variable="intU" tolerance="1e-5">7.5</value>
            <value variable="intV" tolerance="1e-5">1.3790930824</value>
            <value variable="intW" tolerance="1e-5">3</value>
        </metric>
    </metrics>
</test>