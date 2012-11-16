<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Channel Flow boundary conditions from file</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Channel_Flow_2D_BC_from_file.xml</parameters>
    <files>
        <file description="Session File">Channel_Flow_2D_BC_from_file.xml</file>
        <file description="Session File">Channel_Flow_2D_BC_from_file_u_1.bc</file>
        <file description="Session File">Channel_Flow_2D_BC_from_file_u_3.bc</file>
        <file description="Session File">Channel_Flow_2D_BC_from_file_v_3.bc</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">1.85224e-13</value>
            <value variable="v" tolerance="1e-12">2.65241e-13</value>
            <value variable="p" tolerance="1e-8">4.8272e-11</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">8.31335e-13</value>
            <value variable="v" tolerance="1e-12">2.54448e-13</value>
            <value variable="p" tolerance="1e-8">1.03323e-10</value>
        </metric>
    </metrics>
</test>


