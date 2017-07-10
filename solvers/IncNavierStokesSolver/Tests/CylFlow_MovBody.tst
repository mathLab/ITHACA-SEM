<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D flexible cylinder flow simulation using "MovingBody" module</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>CylFlow_MovBody.xml</parameters>
    <files>
        <file description="Session File">CylFlow_MovBody.xml</file>
    </files>
    <metrics>
L 2 error (variable u) : 55.4979
L inf error (variable u) : 1.44442
L 2 error (variable v) : 3.73646
L inf error (variable v) : 0.670941
L 2 error (variable w) : 0.00651794
L inf error (variable w) : 0.00445146
L 2 error (variable p) : 167.909
L inf error (variable p) : 5.6965
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">55.4979</value>
            <value variable="v" tolerance="1e-12">3.73646</value>
            <value variable="w" tolerance="1e-12">0.00651794</value>
            <value variable="p" tolerance="1e-12">167.909</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.44442</value>
            <value variable="v" tolerance="1e-12">0.670941</value>
            <value variable="w" tolerance="1e-12">0.00445146</value>
            <value variable="p" tolerance="1e-12">5.6965</value>
        </metric>
    </metrics>
</test>
