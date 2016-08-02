<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Channel Flow P=8 BodyForce</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>--use-metis ChanFlow_m8_BodyForce.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File">ChanFlow_m8_BodyForce.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
L 2 error (variable u) : 9.90752e-06
L inf error (variable u) : 1.41153e-05
L 2 error (variable v) : 4.56171e-10
L inf error (variable v) : 1.45252e-09
L 2 error (variable p) : 1.50907e-07
L inf error (variable p) : 5.30544e-07
            <value variable="u" tolerance="1e-8">9.90752e-06</value>
            <value variable="v" tolerance="1e-8">4.56171e-10</value>
            <value variable="p" tolerance="1e-8">1.50907e-07</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-8">1.41153e-05</value>
            <value variable="v" tolerance="1e-8">1.45252e-09</value>
            <value variable="p" tolerance="1e-8">5.30544e-07</value>
        </metric>
    </metrics>
</test>
