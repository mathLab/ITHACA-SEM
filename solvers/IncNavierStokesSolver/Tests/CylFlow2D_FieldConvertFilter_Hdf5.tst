<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>2D cylinder flow, P=4, parallel HDF5 input/output and FieldConvert filter</description>    
    <executable>IncNavierStokesSolver</executable>
    <parameters>--io-format Hdf5 CylFlow2D_FieldConvertFilter_Hdf5.xml</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">CylFlow2D_FieldConvertFilter_Hdf5.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-8">21.4085</value>
            <value variable="v" tolerance="1e-8">0.640317</value>
            <value variable="p" tolerance="1e-6">1.03461</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-8">1.84427</value>
            <value variable="v" tolerance="1e-8">0.872692</value>
            <value variable="p" tolerance="1e-6">1.38755</value>
        </metric>
    </metrics>
</test>
