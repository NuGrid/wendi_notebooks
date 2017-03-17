import fileupload
from IPython.display import display
# Install Javascript
fileupload.nbinstall()
def _upload():
    import io
    _upload_widget = fileupload.FileUploadWidget()

    def _cb(change):
        decoded = io.StringIO(change['owner'].data.decode('utf-8'))
        filename = change['owner'].filename
        print('Uploaded `{}` ({:.2f} kB)'.format(
            filename, len(decoded.read()) / 2 **10))
        output= decoded.getvalue()
        filename_out='star_data.txt'
        f1=open(filename_out,'w')
        f1.write(output)
        f1.close()
    _upload_widget.observe(_cb, names='data')
    display(_upload_widget)
