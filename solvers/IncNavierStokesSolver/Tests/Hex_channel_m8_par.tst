<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>3D channel flow, Hex elements, par(2), P=8</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Hex_channel_m8_par.xml</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">Hex_channel_m8_par.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">3.20512e-13</value>
            <value variable="v" tolerance="1e-12">2.4627e-13</value>
            <value variable="w" tolerance="1e-08">3.0802e-12</value>
            <value variable="p" tolerance="1e-08">2.60141e-11</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.17002e-12</value>
            <value variable="v" tolerance="1e-12">8.49471e-13</value>
            <value variable="w" tolerance="1e-08">1.74388e-11</value>
            <value variable="p" tolerance="1e-08">9.24516e-11</value>
        </metric>
    </metrics>
</test>
