<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Convert compressed xml mesh to vtu with scaling </description>
    <executable>FieldConvert</executable>
    <parameters> -f -e scale_compressed_2D.xml scale_compresed_2D.vtu </parameters>
    <files>
        <file description="Mesh File">scale_compressed_2D.xml</file>
    </files>
     <metrics>
        <metric type="L2" id="1">
            <value variable="x"       tolerance="1e-6">3.04817</value>
            <value variable="y"       tolerance="1e-6">2.7107</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="x"       tolerance="1e-6">1.91443</value>
            <value variable="y"       tolerance="1e-6">2.04264</value>
        </metric>
    </metrics>
</test>