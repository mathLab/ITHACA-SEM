<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D cylinder flow, mixed elements, P=5</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Cyl_AdaptiveSFD.xml.gz Cyl_AdaptiveSFD.xml</parameters>
    <files>
        <file description="Session File">Cyl_AdaptiveSFD.xml.gz</file>
        <file description="Session File">Cyl_AdaptiveSFD.xml</file>
        <file description="Session File">Cyl_AdaptiveSFD_LinNS.xml</file>
        <file description="Session File">Cyl_AdaptiveSFD.rst</file>
        <file description="Session File">Cyl_AdaptiveSFD_Eig.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-8">0.241485</value>
            <value variable="v" tolerance="1e-8">0.254726</value>
            <value variable="p" tolerance="1e-4">0.0726467</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-8">0.0829632</value>
            <value variable="v" tolerance="1e-8">0.0752455</value>
            <value variable="p" tolerance="1e-8">0.0251153</value>
        </metric>
    </metrics>
</test>
