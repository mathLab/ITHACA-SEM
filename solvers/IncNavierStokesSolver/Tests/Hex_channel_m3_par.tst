<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D channel flow, Hex elements, par(2), P=3</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>--use-scotch Hex_channel_m3_par.xml</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">Hex_channel_m3_par.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-08">1.03412e-11</value>
            <value variable="v" tolerance="1e-08">7.62679e-12</value>
            <value variable="w" tolerance="1e-08">8.14841e-11</value>
            <value variable="p" tolerance="1e-08">2.21238e-09</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-08">9.31369e-11</value>
            <value variable="v" tolerance="1e-08">7.31505e-11</value>
            <value variable="w" tolerance="1e-08">5.83095e-09</value>
            <value variable="p" tolerance="1e-07">9.42423e-08</value>
        </metric>
    </metrics>
</test>
