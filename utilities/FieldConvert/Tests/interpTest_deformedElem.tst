<test>
    <description>Test interpolation in deformed elems on surfaces</description>
    <executable>FieldConvert</executable>
    <parameters>-m interpfield:fromxml=interpTest_mesh_1.xml:fromfld=interpTest_field_1.fld interpTest_mesh_2.xml interpTest_field_2.fld</parameters>
    <files>
        <file description="Session File">interpTest_mesh_2.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="0.046">0</value>
        </metric>
    </metrics>
</test>

